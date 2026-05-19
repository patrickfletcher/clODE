from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np

from ..observers._definitions import ResolvedObserverSpec, resolve_observer_spec
from ..observers.types import (
    ObserverParams,
    _EventOutputSettings,
)
from ..problem._core import ProblemInfo, RhsSource
from ..simulation.params import (
    SolverParams,
    _TrajectoryOutputSettings,
)
from ..simulation._stepper_definitions import StepperDefinition
from .buffers import (
    ArrayLayout,
    BufferManager,
    CommonBuffers,
    FeatureBuffers,
    N_RNGSTATE,
    TrajectoryBuffers,
)
from .models import KernelKind, Precision, ProblemShape, ProgramBundle
from .observer_metadata import ObserverMetadata, get_observer_metadata
from .registry import KernelRegistry
from .runtime import OpenCLRuntime, _require_opencl_binding
from .source_builder import SourceBuilder


@dataclass(slots=True)
class _TransientTransferCache:
    initial_state: np.ndarray | None = None
    current_dt: np.ndarray | None = None
    final_state: np.ndarray | None = None
    final_time: np.ndarray | None = None
    has_result: bool = False

    def reset_solver_state(
        self,
        *,
        initial_state: np.ndarray | None = None,
        current_dt: np.ndarray | None = None,
    ) -> None:
        if initial_state is not None:
            self.initial_state = np.array(initial_state, dtype=np.float64, copy=True)
        if current_dt is not None:
            self.current_dt = np.array(current_dt, dtype=np.float64, copy=True)
        self.final_state = None
        self.final_time = None
        self.has_result = False

    def invalidate_result(self) -> None:
        self.final_state = None
        self.final_time = None
        self.has_result = False

    def mark_result_pending(self) -> None:
        self.final_state = None
        self.current_dt = None
        self.final_time = None
        self.has_result = True

    def promote_final_state_to_initial_state(self) -> None:
        self.initial_state = (
            None
            if self.final_state is None
            else np.array(self.final_state, dtype=np.float64, copy=True)
        )


@dataclass(slots=True)
class _TrajectoryTransferCache:
    t: np.ndarray | None = None
    x: np.ndarray | None = None
    dx: np.ndarray | None = None
    aux: np.ndarray | None = None
    n_stored: np.ndarray | None = None
    has_result: bool = False

    def invalidate(self) -> None:
        self.t = None
        self.x = None
        self.dx = None
        self.aux = None
        self.n_stored = None
        self.has_result = False

    def mark_result_pending(self) -> None:
        self.t = None
        self.x = None
        self.dx = None
        self.aux = None
        self.n_stored = None
        self.has_result = True


@dataclass(slots=True)
class _FeatureTransferCache:
    features: np.ndarray | None = None
    has_result: bool = False

    def invalidate(self) -> None:
        self.features = None
        self.has_result = False

    def mark_result_pending(self) -> None:
        self.features = None
        self.has_result = True


class OpenCLTransientExecutor:
    def __init__(
        self,
        problem_info: ProblemInfo,
        rhs_source: RhsSource,
        stepper: str,
        single_precision: bool,
        runtime: OpenCLRuntime,
        clode_root: str,
    ) -> None:
        self._problem_info = problem_info
        self._rhs_source = rhs_source
        self._stepper = stepper
        self._precision = Precision.from_single_precision(single_precision)
        self._runtime = runtime
        self._kernel_root = Path(clode_root)
        self._registry = KernelRegistry(self._kernel_root)
        self._stepper_definition: StepperDefinition = self._registry.get_stepper_definition(
            stepper
        )
        self._stepper = self._stepper_definition.stepper_name
        self._source_builder = SourceBuilder(self._kernel_root, self._registry)
        self._buffer_manager = BufferManager(self._runtime, self._precision)
        self._problem_shape = ProblemShape.from_problem_info(problem_info)
        self._opencl_binding = _require_opencl_binding()

        self._solver_params = SolverParams()
        self._integration_settings = self._solver_params.integration_settings
        self._tspan: tuple[float, float] = (0.0, 0.0)
        self._buffers: CommonBuffers | None = None
        self._program_bundle: ProgramBundle | None = None
        self._transient_cache = _TransientTransferCache()
        self._pending_seed: int | None = None

    @staticmethod
    def _copy_solver_params(solver_params: SolverParams) -> SolverParams:
        return solver_params.copy()

    def build_cl(self) -> None:
        if self._precision is Precision.DOUBLE:
            self._runtime.require_double_precision()

        source_bundle = self._source_builder.build(
            kernel_kind=KernelKind.TRANSIENT,
            precision=self._precision,
            stepper_name=self._stepper_definition.stepper_name,
            problem_shape=self._problem_shape,
            rhs=self._rhs_source,
        )
        self._program_bundle = self._runtime.program_cache.get_or_build(
            self._runtime, source_bundle
        )

    def get_available_steppers(self) -> list[str]:
        return list(self._registry.get_available_steppers())

    def get_dt(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._transient_cache.current_dt is None:
            self._transient_cache.current_dt = self._buffer_manager.download_dt(
                self._buffers
            ).reshape(-1)
        return self._transient_cache.current_dt.astype(np.float64, copy=False).tolist()

    def get_program_string(self) -> str:
        if self._program_bundle is None:
            return ""
        source_bundle = self._program_bundle.source_bundle
        options = "\n".join(source_bundle.build_options)
        if not options:
            return source_bundle.source_text
        return f"// Build options:\n// {options.replace(chr(10), chr(10) + '// ')}\n{source_bundle.source_text}"

    def get_solver_params(self) -> SolverParams:
        return self._copy_solver_params(self._solver_params)

    def get_tf(self) -> list[float]:
        if self._buffers is None or not self._transient_cache.has_result:
            return []
        if self._transient_cache.final_time is None:
            self._transient_cache.final_time = self._buffer_manager.download_tf(
                self._buffers
            ).reshape(-1)
        return self._transient_cache.final_time.astype(
            np.float64, copy=False
        ).tolist()

    def get_tspan(self) -> list[float]:
        return [float(self._tspan[0]), float(self._tspan[1])]

    def get_x0(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._transient_cache.initial_state is None:
            x0 = self._buffer_manager.download_x0(self._buffers)
            self._transient_cache.initial_state = x0.flatten(order="F")
        return self._transient_cache.initial_state.astype(
            np.float64, copy=False
        ).tolist()

    def get_xf(self) -> list[float] | None:
        if self._buffers is None or not self._transient_cache.has_result:
            return None
        if self._transient_cache.final_state is None:
            xf = self._buffer_manager.download_xf(self._buffers)
            self._transient_cache.final_state = xf.flatten(order="F")
        return self._transient_cache.final_state.astype(
            np.float64, copy=False
        ).tolist()

    def print_status(self) -> None:
        print("------------------")
        print(f"   {self._rhs_source.origin_label}")
        print(f"   nVar={self._problem_shape.n_var}")
        print(f"   nPar={self._problem_shape.n_par}")
        print(f"   nAux={self._problem_shape.n_aux}")
        print(f"   nWiener={self._problem_shape.n_wiener}")
        print(
            f"Using {'single' if self._precision is Precision.SINGLE else 'double'} precision."
        )
        print(f"Using stepper: {self._stepper}")
        print(f"Using OpenCL runtime: {self._runtime.describe()}")

    def seed_rng(self, seed: int | None = None) -> None:
        if self._buffers is None:
            self._pending_seed = seed
            return

        rng_state = self._make_rng_state(self._buffers.ensemble_size, seed)
        self._buffer_manager.upload_rng_state(self._buffers, rng_state)
        self._buffer_manager.reset_rng_box_muller_cache(self._buffers)
        self._buffer_manager.clear_prepared_wiener_state(self._buffers)

    def set_pars(self, parameters: Sequence[float]) -> None:
        buffers = self._require_buffers()
        expected = buffers.ensemble_size * self._problem_shape.n_par
        host = np.asarray(parameters, dtype=np.float64)
        if host.size != expected:
            raise ValueError(
                f"Expected {expected} parameter values, received {host.size}"
            )
        matrix = host.reshape(
            (buffers.ensemble_size, self._problem_shape.n_par), order="F"
        )
        self._buffer_manager.upload_pars(buffers, matrix)
        self._reset_solver_state_cache()

    def set_problem_data(
        self, initial_state: Sequence[float], parameters: Sequence[float]
    ) -> None:
        initial_state_host = np.asarray(initial_state, dtype=np.float64)
        parameter_host = np.asarray(parameters, dtype=np.float64)
        ensemble_size = self._infer_ensemble_size(initial_state_host, parameter_host)

        buffers_reallocated = (
            self._buffers is None or self._buffers.ensemble_size != ensemble_size
        )
        if buffers_reallocated:
            self._buffers = self._buffer_manager.allocate_common(
                ensemble_size=ensemble_size,
                shape=self._problem_shape,
            )
            self._buffer_manager.upload_tspan(self._buffers, self._tspan)
            self._buffer_manager.upload_integration_settings(
                self._buffers,
                self._integration_settings,
            )

        x0_matrix = initial_state_host.reshape(
            (ensemble_size, self._problem_shape.n_var), order="F"
        )
        if self._problem_shape.n_par > 0:
            pars_matrix = parameter_host.reshape(
                (ensemble_size, self._problem_shape.n_par), order="F"
            )
        else:
            pars_matrix = np.empty((ensemble_size, 0), dtype=np.float64)

        self._buffer_manager.upload_problem_data(self._buffers, x0_matrix, pars_matrix)
        self._buffer_manager.clear_prepared_wiener_state(self._buffers)
        self._reset_solver_state_cache(initial_state=initial_state_host)

        if buffers_reallocated:
            pending_seed = self._pending_seed
            self._pending_seed = None
            self.seed_rng(pending_seed)

    def set_solver_params(self, solver_params: SolverParams) -> None:
        self._solver_params = self._copy_solver_params(solver_params)
        self._integration_settings = self._solver_params.integration_settings
        self._apply_integration_settings_to_runtime(reset_solver_state=True)

    def set_tspan(self, tspan: Sequence[float]) -> None:
        if len(tspan) != 2:
            raise ValueError("tspan must contain exactly two values")
        self._tspan = (float(tspan[0]), float(tspan[1]))
        if self._buffers is not None:
            self._buffer_manager.upload_tspan(self._buffers, self._tspan)
        self._transient_cache.invalidate_result()

    def set_x0(self, initial_state: Sequence[float]) -> None:
        buffers = self._require_buffers()
        expected = buffers.ensemble_size * self._problem_shape.n_var
        host = np.asarray(initial_state, dtype=np.float64)
        if host.size != expected:
            raise ValueError(
                f"Expected {expected} initial-state values, received {host.size}"
            )
        matrix = host.reshape(
            (buffers.ensemble_size, self._problem_shape.n_var), order="F"
        )
        self._buffer_manager.upload_x0(buffers, matrix)
        self._buffer_manager.clear_prepared_wiener_state(buffers)
        self._reset_solver_state_cache(initial_state=host)

    def shift_tspan(self) -> None:
        start, end = self._tspan
        duration = end - start
        self.set_tspan((end, end + duration))

    def shift_x0(self) -> None:
        buffers = self._require_buffers()
        self._opencl_binding.enqueue_copy(
            self._runtime.queue,
            buffers.x0,
            buffers.xf,
            byte_count=buffers.x0.size,
            src_offset=0,
            dst_offset=0,
        ).wait()
        self._transient_cache.promote_final_state_to_initial_state()

    def transient(self) -> None:
        if self._program_bundle is None:
            raise RuntimeError("OpenCL program has not been built")
        buffers = self._require_buffers()
        kernel = self._program_bundle.kernels["transient"]
        kernel.set_args(
            buffers.tspan,
            buffers.x0,
            buffers.pars,
            buffers.integration_settings,
            buffers.xf,
            buffers.rng_state,
            buffers.rng_spare_normal,
            buffers.rng_spare_normal_valid,
            buffers.prepared_wiener,
            buffers.prepared_wiener_valid,
            buffers.dt,
            buffers.tf,
        )
        self._opencl_binding.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._transient_cache.mark_result_pending()

    def _reset_solver_state_cache(
        self,
        *,
        initial_state: np.ndarray | None = None,
    ) -> None:
        buffers = self._require_buffers()
        current_dt = self._buffer_manager.reset_solver_dt(
            buffers,
            self._integration_settings,
        )
        self._transient_cache.reset_solver_state(
            initial_state=initial_state,
            current_dt=current_dt.reshape(-1),
        )

    def _apply_integration_settings_to_runtime(
        self, *, reset_solver_state: bool
    ) -> None:
        if self._buffers is None:
            return

        self._buffer_manager.upload_integration_settings(
            self._buffers,
            self._integration_settings,
        )
        if reset_solver_state:
            self._buffer_manager.clear_prepared_wiener_state(self._buffers)
            self._reset_solver_state_cache()

    def _infer_ensemble_size(
        self, initial_state: np.ndarray, parameters: np.ndarray
    ) -> int:
        if self._problem_shape.n_var == 0:
            raise ValueError("Problem must contain at least one variable")
        if initial_state.size % self._problem_shape.n_var != 0:
            raise ValueError(
                "Initial-state vector length must be a multiple of the variable count"
            )

        ensemble_size = initial_state.size // self._problem_shape.n_var
        if self._problem_shape.n_par == 0:
            if parameters.size != 0:
                raise ValueError(
                    "Parameter vector must be empty when the problem has no parameters"
                )
            return ensemble_size

        if parameters.size % self._problem_shape.n_par != 0:
            raise ValueError(
                "Parameter vector length must be a multiple of the parameter count"
            )
        parameter_ensemble_size = parameters.size // self._problem_shape.n_par
        if parameter_ensemble_size != ensemble_size:
            raise ValueError(
                "Initial-state and parameter vectors must describe the same ensemble size"
            )
        return ensemble_size

    def _make_rng_state(self, ensemble_size: int, seed: int | None) -> np.ndarray:
        if seed is None:
            generator = np.random.default_rng()
            return generator.integers(
                0,
                np.iinfo(np.uint64).max,
                size=(ensemble_size, N_RNGSTATE),
                dtype=np.uint64,
            )

        values = np.arange(
            seed,
            seed + ensemble_size * N_RNGSTATE,
            dtype=np.uint64,
        )
        return values.reshape((ensemble_size, N_RNGSTATE), order="F")

    def _require_buffers(self) -> CommonBuffers:
        if self._buffers is None:
            raise RuntimeError("Problem data has not been initialized")
        return self._buffers


class OpenCLTrajectoryExecutor(OpenCLTransientExecutor):
    def __init__(
        self,
        problem_info: ProblemInfo,
        rhs_source: RhsSource,
        stepper: str,
        single_precision: bool,
        runtime: OpenCLRuntime,
        clode_root: str,
    ) -> None:
        super().__init__(
            problem_info,
            rhs_source,
            stepper,
            single_precision,
            runtime,
            clode_root,
        )
        self._trajectory_buffers: TrajectoryBuffers | None = None
        self._trajectory_cache = _TrajectoryTransferCache()
        self._trajectory_output_settings = self._solver_params.trajectory_output_settings

    def build_cl(self) -> None:
        if self._precision is Precision.DOUBLE:
            self._runtime.require_double_precision()

        source_bundle = self._source_builder.build(
            kernel_kind=KernelKind.TRAJECTORY,
            precision=self._precision,
            stepper_name=self._stepper_definition.stepper_name,
            problem_shape=self._problem_shape,
            rhs=self._rhs_source,
        )
        self._program_bundle = self._runtime.program_cache.get_or_build(
            self._runtime, source_bundle
        )

    def set_problem_data(
        self, initial_state: Sequence[float], parameters: Sequence[float]
    ) -> None:
        previous_ensemble_size = None if self._buffers is None else self._buffers.ensemble_size
        super().set_problem_data(initial_state, parameters)
        current_ensemble_size = self._require_buffers().ensemble_size
        if previous_ensemble_size != current_ensemble_size and self._trajectory_buffers is not None:
            self._trajectory_buffers = None

    def set_solver_params(self, solver_params: SolverParams) -> None:
        previous_integration = self._integration_settings
        previous_output = self._trajectory_output_settings

        self._solver_params = self._copy_solver_params(solver_params)
        self._integration_settings = self._solver_params.integration_settings
        self._trajectory_output_settings = self._solver_params.trajectory_output_settings

        integration_changed = previous_integration != self._integration_settings
        output_changed = previous_output != self._trajectory_output_settings

        if not integration_changed and not output_changed:
            return

        self._apply_integration_settings_to_runtime(
            reset_solver_state=integration_changed
        )

        if output_changed and self._trajectory_buffers is not None:
            if previous_output.max_store != self._trajectory_output_settings.max_store:
                self._trajectory_buffers = None
            else:
                self._buffer_manager.upload_trajectory_output_settings(
                    self._trajectory_buffers,
                    self._trajectory_output_settings,
                )

        if integration_changed or output_changed:
            self._invalidate_trajectory_cache()

    def get_aux(self) -> list[float]:
        if self._buffers is None or not self._trajectory_cache.has_result:
            return []
        if self._trajectory_cache.aux is None:
            self._trajectory_cache.aux = self._buffer_manager.download_aux(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._trajectory_cache.aux.astype(np.float64, copy=False).tolist()

    def get_dx(self) -> list[float]:
        if self._buffers is None or not self._trajectory_cache.has_result:
            return []
        if self._trajectory_cache.dx is None:
            self._trajectory_cache.dx = self._buffer_manager.download_dx(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._trajectory_cache.dx.astype(np.float64, copy=False).tolist()

    def get_n_stored(self) -> list[int]:
        if self._buffers is None or not self._trajectory_cache.has_result:
            return []
        if self._trajectory_cache.n_stored is None:
            self._trajectory_cache.n_stored = self._buffer_manager.download_n_stored(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._trajectory_cache.n_stored.astype(np.int32, copy=False).tolist()

    def get_t(self) -> list[float]:
        if self._buffers is None or not self._trajectory_cache.has_result:
            return []
        if self._trajectory_cache.t is None:
            self._trajectory_cache.t = self._buffer_manager.download_t(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._trajectory_cache.t.astype(np.float64, copy=False).tolist()

    def get_x(self) -> list[float]:
        if self._buffers is None or not self._trajectory_cache.has_result:
            return []
        if self._trajectory_cache.x is None:
            self._trajectory_cache.x = self._buffer_manager.download_x(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._trajectory_cache.x.astype(np.float64, copy=False).tolist()

    def trajectory(self) -> None:
        if self._program_bundle is None:
            raise RuntimeError("OpenCL program has not been built")
        buffers = self._require_buffers()
        trajectory_buffers = self._ensure_trajectory_buffers()
        kernel = self._program_bundle.kernels["trajectory"]
        kernel.set_args(
            buffers.tspan,
            buffers.x0,
            buffers.pars,
            buffers.integration_settings,
            buffers.xf,
            buffers.rng_state,
            buffers.rng_spare_normal,
            buffers.rng_spare_normal_valid,
            buffers.prepared_wiener,
            buffers.prepared_wiener_valid,
            buffers.dt,
            buffers.tf,
            trajectory_buffers.output_settings,
            trajectory_buffers.t,
            trajectory_buffers.x,
            trajectory_buffers.dx,
            trajectory_buffers.aux,
            trajectory_buffers.n_stored,
        )
        self._opencl_binding.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._transient_cache.mark_result_pending()
        self._trajectory_cache.mark_result_pending()

    def _ensure_trajectory_buffers(self) -> TrajectoryBuffers:
        buffers = self._require_buffers()
        if self._trajectory_buffers is None:
            self._trajectory_buffers = self._buffer_manager.allocate_trajectory(
                buffers.ensemble_size,
                buffers.problem_shape,
                self._trajectory_output_settings,
            )
            self._buffer_manager.upload_trajectory_output_settings(
                self._trajectory_buffers,
                self._trajectory_output_settings,
            )
        return self._trajectory_buffers

    def _invalidate_trajectory_cache(self) -> None:
        self._trajectory_cache.invalidate()

    def _require_trajectory_buffers(self) -> TrajectoryBuffers:
        if self._trajectory_buffers is None:
            raise RuntimeError("Trajectory buffers have not been initialized")
        return self._trajectory_buffers


class OpenCLFeatureExecutor(OpenCLTransientExecutor):
    def __init__(
        self,
        problem_info: ProblemInfo,
        rhs_source: RhsSource,
        stepper: str,
        observer: str,
        observer_params: ObserverParams,
        single_precision: bool,
        runtime: OpenCLRuntime,
        clode_root: str,
    ) -> None:
        super().__init__(
            problem_info,
            rhs_source,
            stepper,
            single_precision,
            runtime,
            clode_root,
        )
        self._observer_name = observer
        self._observer_params = self._copy_observer_params(observer_params)
        self._observer_runtime_settings = observer_params.runtime_settings
        self._event_output_settings = observer_params.event_output_settings
        self._resolved_observer_spec = self._resolve_observer_spec()
        self._feature_metadata = self._resolve_feature_metadata()
        self._feature_buffers: FeatureBuffers | None = None
        self._feature_cache = _FeatureTransferCache()
        self._observer_initialized = False

    @staticmethod
    def _copy_observer_params(observer_params: ObserverParams) -> ObserverParams:
        return ObserverParams(
            observer_params.e_var_ix,
            observer_params.f_var_ix,
            observer_params.max_event_count,
            observer_params.max_event_timestamps,
            observer_params.min_amp,
            observer_params.min_imi,
            observer_params.nhood_radius,
            observer_params.x_up_threshold,
            observer_params.x_down_threshold,
            observer_params.dx_up_threshold,
            observer_params.dx_down_threshold,
            observer_params.eps_dx,
        )

    def build_cl(self) -> None:
        if self._precision is Precision.DOUBLE:
            self._runtime.require_double_precision()

        source_bundle = self._source_builder.build(
            kernel_kind=KernelKind.FEATURES,
            precision=self._precision,
            stepper_name=self._stepper_definition.stepper_name,
            problem_shape=self._problem_shape,
            rhs=self._rhs_source,
            resolved_observer_spec=self._resolved_observer_spec,
        )
        self._program_bundle = self._runtime.program_cache.get_or_build(
            self._runtime, source_bundle
        )

    def features(self, reinitialize_observer: bool | None = None) -> None:
        if self._program_bundle is None:
            raise RuntimeError("OpenCL program has not been built")
        if reinitialize_observer is True or not self._observer_initialized:
            self.initialize_observer()

        buffers = self._require_buffers()
        feature_buffers = self._ensure_feature_buffers()
        kernel = self._program_bundle.kernels["features"]
        kernel.set_args(
            buffers.tspan,
            buffers.x0,
            buffers.pars,
            buffers.integration_settings,
            buffers.xf,
            buffers.rng_state,
            buffers.rng_spare_normal,
            buffers.rng_spare_normal_valid,
            buffers.prepared_wiener,
            buffers.prepared_wiener_valid,
            buffers.dt,
            buffers.tf,
            feature_buffers.observer_state,
            feature_buffers.observer_runtime_settings,
            feature_buffers.features,
        )
        self._opencl_binding.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._transient_cache.mark_result_pending()
        self._observer_initialized = True
        self._feature_cache.mark_result_pending()

    def get_f(self) -> list[float]:
        if self._buffers is None or not self._feature_cache.has_result:
            return []
        if self._feature_cache.features is None:
            self._feature_cache.features = self._buffer_manager.download_features(
                self._buffers,
                self._require_feature_buffers(),
            )
        return self._feature_cache.features.astype(np.float64, copy=False).tolist()

    def get_feature_names(self) -> list[str]:
        return list(self._feature_metadata.feature_names)

    def get_n_features(self) -> int:
        return self._feature_metadata.n_features

    def get_observer_params(self) -> ObserverParams:
        return self._copy_observer_params(self._observer_params)

    def initialize_observer(self) -> None:
        if self._program_bundle is None:
            raise RuntimeError("OpenCL program has not been built")
        buffers = self._require_buffers()
        feature_buffers = self._ensure_feature_buffers()
        self._buffer_manager.upload_observer_runtime_settings(
            feature_buffers,
            self._observer_runtime_settings,
        )
        self._buffer_manager.clear_observer_state(feature_buffers, buffers.ensemble_size)
        kernel = self._program_bundle.kernels["initializeObserver"]
        kernel.set_args(
            buffers.tspan,
            buffers.x0,
            buffers.pars,
            buffers.integration_settings,
            buffers.rng_state,
            buffers.rng_spare_normal,
            buffers.rng_spare_normal_valid,
            buffers.prepared_wiener,
            buffers.prepared_wiener_valid,
            buffers.dt,
            feature_buffers.observer_state,
            feature_buffers.observer_runtime_settings,
        )
        self._opencl_binding.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._observer_initialized = True
        self._feature_cache.invalidate()

    def is_observer_initialized(self) -> bool:
        return self._observer_initialized

    def set_observer(self, observer: str) -> None:
        if observer == self._observer_name:
            return
        self._registry.validate_observer(observer)
        self._observer_name = observer
        self._program_bundle = None
        self._observer_initialized = False
        self._invalidate_feature_cache()
        self._recreate_feature_metadata()

    def set_observer_params(self, observer_params: ObserverParams) -> None:
        previous_spec = self._resolved_observer_spec
        previous_runtime_settings = self._observer_runtime_settings
        previous_metadata = self._feature_metadata
        self._observer_params = self._copy_observer_params(observer_params)
        self._observer_runtime_settings = self._observer_params.runtime_settings
        self._event_output_settings = _EventOutputSettings(
            self._observer_params.max_event_timestamps
        )
        self._resolved_observer_spec = self._resolve_observer_spec()
        build_changed = (
            previous_spec.build_signature != self._resolved_observer_spec.build_signature
        )
        if build_changed:
            self._program_bundle = None
        self._observer_initialized = False
        self._invalidate_feature_cache()
        self._feature_metadata = self._resolve_feature_metadata()

        if self._feature_buffers is not None:
            buffer_shape_changed = (
                self._resolved_observer_spec.layout_signature
                != previous_spec.layout_signature
                or self._feature_metadata.n_features != previous_metadata.n_features
            )
            if buffer_shape_changed:
                self._feature_buffers = None
            elif previous_runtime_settings != self._observer_runtime_settings:
                self._buffer_manager.upload_observer_runtime_settings(
                    self._feature_buffers,
                    self._observer_runtime_settings,
                )

    def set_pars(self, parameters: Sequence[float]) -> None:
        super().set_pars(parameters)
        self._observer_initialized = False
        self._invalidate_feature_cache()

    def set_problem_data(
        self, initial_state: Sequence[float], parameters: Sequence[float]
    ) -> None:
        previous_ensemble_size = None if self._buffers is None else self._buffers.ensemble_size
        super().set_problem_data(initial_state, parameters)
        current_ensemble_size = self._require_buffers().ensemble_size
        if previous_ensemble_size != current_ensemble_size and self._feature_buffers is not None:
            self._feature_buffers = None
        self._observer_initialized = False
        self._invalidate_feature_cache()

    def set_solver_params(self, solver_params: SolverParams) -> None:
        super().set_solver_params(solver_params)
        self._invalidate_feature_cache()

    def set_tspan(self, tspan: Sequence[float]) -> None:
        super().set_tspan(tspan)
        self._invalidate_feature_cache()

    def set_x0(self, initial_state: Sequence[float]) -> None:
        super().set_x0(initial_state)
        self._observer_initialized = False
        self._invalidate_feature_cache()

    def _ensure_feature_buffers(self) -> FeatureBuffers:
        buffers = self._require_buffers()
        if self._feature_buffers is None:
            self._feature_buffers = self._buffer_manager.allocate_feature(
                buffers.ensemble_size,
                self._feature_metadata.n_features,
                self._feature_metadata.observer_state_nbytes,
            )
            self._buffer_manager.upload_observer_runtime_settings(
                self._feature_buffers,
                self._observer_runtime_settings,
            )
            self._buffer_manager.clear_observer_state(
                self._feature_buffers,
                buffers.ensemble_size,
            )
        return self._feature_buffers

    def _invalidate_feature_cache(self) -> None:
        self._feature_cache.invalidate()

    def _recreate_feature_metadata(self) -> None:
        if self._feature_buffers is not None:
            self._feature_buffers = None
        self._resolved_observer_spec = self._resolve_observer_spec()
        self._feature_metadata = self._resolve_feature_metadata()

    def _require_feature_buffers(self) -> FeatureBuffers:
        if self._feature_buffers is None:
            raise RuntimeError("Feature buffers have not been initialized")
        return self._feature_buffers

    def _resolve_feature_metadata(self) -> ObserverMetadata:
        return get_observer_metadata(
            self._runtime,
            self._problem_info,
            self._observer_name,
            self._observer_params,
            self._precision,
            resolved_observer_spec=self._resolved_observer_spec,
        )

    def _resolve_observer_spec(self) -> ResolvedObserverSpec:
        return resolve_observer_spec(
            self._problem_info,
            self._observer_name,
            self._observer_params,
            real_dtype=np.dtype(
                np.float32 if self._precision is Precision.SINGLE else np.float64
            ),
            n_store_events=self._event_output_settings.max_event_timestamps,
        )
