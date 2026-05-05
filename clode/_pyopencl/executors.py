from __future__ import annotations

from pathlib import Path
from typing import Sequence

import numpy as np

from .._backends.protocol import FeatureBackend, SimulatorBackend
from .._backends.rhs import RhsSource
from ..types import ObserverParams, ProblemInfo, SolverParams
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
from .runtime import OpenCLRuntime, _require_pyopencl
from .source_builder import SourceBuilder


class PyOpenCLTransientBackend(SimulatorBackend):
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
        self._source_builder = SourceBuilder(self._kernel_root, self._registry)
        self._buffer_manager = BufferManager(self._runtime, self._precision)
        self._problem_shape = ProblemShape.from_problem_info(problem_info)
        self._pyopencl = _require_pyopencl()

        self._solver_params = SolverParams()
        self._tspan: tuple[float, float] = (0.0, 0.0)
        self._buffers: CommonBuffers | None = None
        self._program_bundle: ProgramBundle | None = None
        self._retired_common_buffers: list[CommonBuffers] = []

        self._x0_host: np.ndarray | None = None
        self._pars_host: np.ndarray | None = None
        self._xf_host: np.ndarray | None = None
        self._dt_host: np.ndarray | None = None
        self._tf_host: np.ndarray | None = None
        self._rng_state_host: np.ndarray | None = None

        self._has_transient_result = False
        self._pending_seed: int | None = None

    @staticmethod
    def _copy_solver_params(solver_params: SolverParams) -> SolverParams:
        return SolverParams(
            solver_params.dt,
            solver_params.dtmax,
            solver_params.abstol,
            solver_params.reltol,
            solver_params.max_steps,
            solver_params.max_store,
            solver_params.nout,
        )

    def build_cl(self) -> None:
        if self._precision is Precision.DOUBLE:
            self._runtime.require_double_precision()

        source_bundle = self._source_builder.build(
            kernel_kind=KernelKind.TRANSIENT,
            precision=self._precision,
            stepper_name=self._stepper,
            problem_shape=self._problem_shape,
            rhs=self._rhs_source,
        )
        self._program_bundle = self._runtime.program_cache.get_or_build(
            self._runtime, source_bundle
        )

    def get_available_steppers(self) -> list[str]:
        return list(self._registry._STEPPERS.keys())

    def get_dt(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._dt_host is None:
            self._dt_host = self._buffer_manager.download_dt(self._buffers).reshape(-1)
        return self._dt_host.astype(np.float64, copy=False).tolist()

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
        if self._buffers is None or not self._has_transient_result:
            return []
        if self._tf_host is None:
            self._tf_host = self._buffer_manager.download_tf(self._buffers).reshape(-1)
        return self._tf_host.astype(np.float64, copy=False).tolist()

    def get_tspan(self) -> list[float]:
        return [float(self._tspan[0]), float(self._tspan[1])]

    def get_x0(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._x0_host is None:
            x0 = self._buffer_manager.download_x0(self._buffers)
            self._x0_host = x0.flatten(order="F")
        return self._x0_host.astype(np.float64, copy=False).tolist()

    def get_xf(self) -> list[float] | None:
        if self._buffers is None or not self._has_transient_result:
            return None
        if self._xf_host is None:
            xf = self._buffer_manager.download_xf(self._buffers)
            self._xf_host = xf.flatten(order="F")
        return self._xf_host.astype(np.float64, copy=False).tolist()

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
        self._rng_state_host = rng_state
        self._buffer_manager.upload_rng_state(self._buffers, rng_state)

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
        self._pars_host = host.copy()
        self._has_transient_result = False

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
            if self._buffers is not None:
                self._retired_common_buffers.append(self._buffers)
            self._buffers = self._buffer_manager.allocate_common(
                ensemble_size=ensemble_size,
                shape=self._problem_shape,
            )
            self._buffer_manager.upload_tspan(self._buffers, self._tspan)
            self._buffer_manager.upload_solver_params(self._buffers, self._solver_params)
            self._dt_host = np.full(ensemble_size, self._solver_params.dt, dtype=np.float64)
            self._buffer_manager.upload_dt(self._buffers, self._dt_host)

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
        self._x0_host = initial_state_host.copy()
        self._pars_host = parameter_host.copy()
        self._xf_host = None
        self._tf_host = None
        self._has_transient_result = False

        if buffers_reallocated:
            pending_seed = self._pending_seed
            self._pending_seed = None
            self.seed_rng(pending_seed)

    def set_solver_params(self, solver_params: SolverParams) -> None:
        self._solver_params = self._copy_solver_params(solver_params)
        if self._buffers is None:
            return

        self._buffer_manager.upload_solver_params(self._buffers, self._solver_params)
        self._dt_host = np.full(
            self._buffers.ensemble_size, self._solver_params.dt, dtype=np.float64
        )
        self._buffer_manager.upload_dt(self._buffers, self._dt_host)
        self._tf_host = None
        self._has_transient_result = False

    def set_tspan(self, tspan: Sequence[float]) -> None:
        if len(tspan) != 2:
            raise ValueError("tspan must contain exactly two values")
        self._tspan = (float(tspan[0]), float(tspan[1]))
        if self._buffers is not None:
            self._buffer_manager.upload_tspan(self._buffers, self._tspan)
        self._tf_host = None
        self._has_transient_result = False

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
        self._x0_host = host.copy()
        self._has_transient_result = False

    def shift_tspan(self) -> None:
        start, end = self._tspan
        duration = end - start
        self.set_tspan((end, end + duration))

    def shift_x0(self) -> None:
        buffers = self._require_buffers()
        self._pyopencl.enqueue_copy(
            self._runtime.queue,
            buffers.x0,
            buffers.xf,
            byte_count=buffers.x0.size,
            src_offset=0,
            dst_offset=0,
        ).wait()
        self._x0_host = None if self._xf_host is None else self._xf_host.copy()

    def transient(self) -> None:
        if self._program_bundle is None:
            raise RuntimeError("OpenCL program has not been built")
        buffers = self._require_buffers()
        kernel = self._program_bundle.kernels["transient"]
        kernel.set_args(
            buffers.tspan,
            buffers.x0,
            buffers.pars,
            buffers.solver_params,
            buffers.xf,
            buffers.rng_state,
            buffers.dt,
            buffers.tf,
        )
        self._pyopencl.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._xf_host = None
        self._dt_host = None
        self._tf_host = None
        self._rng_state_host = None
        self._has_transient_result = True

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


class PyOpenCLTrajectoryBackend(PyOpenCLTransientBackend):
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
        self._retired_trajectory_buffers: list[TrajectoryBuffers] = []
        self._t_host: np.ndarray | None = None
        self._x_host: np.ndarray | None = None
        self._dx_host: np.ndarray | None = None
        self._aux_host: np.ndarray | None = None
        self._n_stored_host: np.ndarray | None = None

    def build_cl(self) -> None:
        if self._precision is Precision.DOUBLE:
            self._runtime.require_double_precision()

        source_bundle = self._source_builder.build(
            kernel_kind=KernelKind.TRAJECTORY,
            precision=self._precision,
            stepper_name=self._stepper,
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
            self._retired_trajectory_buffers.append(self._trajectory_buffers)
            self._trajectory_buffers = None

    def set_solver_params(self, solver_params: SolverParams) -> None:
        previous_max_store = self._solver_params.max_store
        super().set_solver_params(solver_params)
        if previous_max_store != self._solver_params.max_store and self._trajectory_buffers is not None:
            self._retired_trajectory_buffers.append(self._trajectory_buffers)
            self._trajectory_buffers = None
        self._invalidate_trajectory_cache()

    def get_aux(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._aux_host is None:
            self._aux_host = self._buffer_manager.download_aux(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._aux_host.astype(np.float64, copy=False).tolist()

    def get_dx(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._dx_host is None:
            self._dx_host = self._buffer_manager.download_dx(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._dx_host.astype(np.float64, copy=False).tolist()

    def get_n_stored(self) -> list[int]:
        if self._buffers is None:
            return []
        if self._n_stored_host is None:
            self._n_stored_host = self._buffer_manager.download_n_stored(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._n_stored_host.astype(np.int32, copy=False).tolist()

    def get_t(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._t_host is None:
            self._t_host = self._buffer_manager.download_t(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._t_host.astype(np.float64, copy=False).tolist()

    def get_x(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._x_host is None:
            self._x_host = self._buffer_manager.download_x(
                self._buffers, self._require_trajectory_buffers()
            )
        return self._x_host.astype(np.float64, copy=False).tolist()

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
            buffers.solver_params,
            buffers.xf,
            buffers.rng_state,
            buffers.dt,
            buffers.tf,
            trajectory_buffers.t,
            trajectory_buffers.x,
            trajectory_buffers.dx,
            trajectory_buffers.aux,
            trajectory_buffers.n_stored,
        )
        self._pyopencl.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._xf_host = None
        self._dt_host = None
        self._tf_host = None
        self._rng_state_host = None
        self._has_transient_result = True
        self._invalidate_trajectory_cache()

    def _ensure_trajectory_buffers(self) -> TrajectoryBuffers:
        buffers = self._require_buffers()
        if self._trajectory_buffers is None:
            self._trajectory_buffers = self._buffer_manager.allocate_trajectory(
                buffers.ensemble_size,
                buffers.problem_shape,
                self._solver_params.max_store,
            )
        return self._trajectory_buffers

    def _invalidate_trajectory_cache(self) -> None:
        self._t_host = None
        self._x_host = None
        self._dx_host = None
        self._aux_host = None
        self._n_stored_host = None

    def _require_trajectory_buffers(self) -> TrajectoryBuffers:
        if self._trajectory_buffers is None:
            raise RuntimeError("Trajectory buffers have not been initialized")
        return self._trajectory_buffers


class PyOpenCLFeatureBackend(PyOpenCLTransientBackend, FeatureBackend):
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
        self._feature_metadata = self._resolve_feature_metadata()
        self._feature_buffers: FeatureBuffers | None = None
        self._retired_feature_buffers: list[FeatureBuffers] = []
        self._features_host: np.ndarray | None = None
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
            stepper_name=self._stepper,
            problem_shape=self._problem_shape,
            rhs=self._rhs_source,
            observer_name=self._observer_name,
            n_store_events=self._observer_params.max_event_timestamps,
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
            buffers.solver_params,
            buffers.xf,
            buffers.rng_state,
            buffers.dt,
            buffers.tf,
            feature_buffers.observer_data,
            feature_buffers.observer_params,
            feature_buffers.features,
        )
        self._pyopencl.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._xf_host = None
        self._dt_host = None
        self._tf_host = None
        self._rng_state_host = None
        self._has_transient_result = True
        self._observer_initialized = True
        self._features_host = None

    def get_f(self) -> list[float]:
        if self._buffers is None:
            return []
        if self._features_host is None:
            self._features_host = self._buffer_manager.download_features(
                self._buffers,
                self._require_feature_buffers(),
            )
        return self._features_host.astype(np.float64, copy=False).tolist()

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
        self._buffer_manager.upload_observer_params(feature_buffers, self._observer_params)
        self._buffer_manager.clear_observer_data(feature_buffers, buffers.ensemble_size)
        kernel = self._program_bundle.kernels["initializeObserver"]
        kernel.set_args(
            buffers.tspan,
            buffers.x0,
            buffers.pars,
            buffers.solver_params,
            buffers.rng_state,
            buffers.dt,
            feature_buffers.observer_data,
            feature_buffers.observer_params,
        )
        self._pyopencl.enqueue_nd_range_kernel(
            self._runtime.queue,
            kernel,
            (buffers.ensemble_size,),
            None,
        ).wait()
        self._runtime.queue.finish()
        self._observer_initialized = True
        self._features_host = None

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
        previous_max_event_timestamps = self._observer_params.max_event_timestamps
        self._observer_params = self._copy_observer_params(observer_params)
        if previous_max_event_timestamps != self._observer_params.max_event_timestamps:
            self._program_bundle = None
        self._observer_initialized = False
        self._invalidate_feature_cache()
        self._recreate_feature_metadata()

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
            self._retired_feature_buffers.append(self._feature_buffers)
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
                self._feature_metadata.observer_data_nbytes,
            )
            self._buffer_manager.upload_observer_params(
                self._feature_buffers,
                self._observer_params,
            )
            self._buffer_manager.clear_observer_data(
                self._feature_buffers,
                buffers.ensemble_size,
            )
        return self._feature_buffers

    def _invalidate_feature_cache(self) -> None:
        self._features_host = None

    def _recreate_feature_metadata(self) -> None:
        if self._feature_buffers is not None:
            self._retired_feature_buffers.append(self._feature_buffers)
            self._feature_buffers = None
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
        )
