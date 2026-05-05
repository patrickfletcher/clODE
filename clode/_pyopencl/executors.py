from __future__ import annotations

from pathlib import Path
from typing import Sequence

from clode.cpp.clode_cpp_wrapper import ProblemInfo, SolverParams
import numpy as np

from .._backends.protocol import SimulatorBackend
from .._backends.rhs import RhsSource
from .buffers import ArrayLayout, BufferManager, CommonBuffers, N_RNGSTATE
from .models import KernelKind, Precision, ProblemShape, ProgramBundle
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

        self._x0_host: np.ndarray | None = None
        self._pars_host: np.ndarray | None = None
        self._xf_host: np.ndarray | None = None
        self._dt_host: np.ndarray | None = None
        self._tf_host: np.ndarray | None = None
        self._rng_state_host: np.ndarray | None = None

        self._has_transient_result = False
        self._pending_seed: int | None = None

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
        return self._program_bundle.source_bundle.source_text

    def get_solver_params(self) -> SolverParams:
        return self._solver_params

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
        self._solver_params = solver_params
        if self._buffers is None:
            return

        self._buffer_manager.upload_solver_params(self._buffers, solver_params)
        self._dt_host = np.full(
            self._buffers.ensemble_size, solver_params.dt, dtype=np.float64
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
