from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..observers.types import ObserverParams
from ..simulation.params import SolverParams
from .models import Precision, ProblemShape
from .runtime import OpenCLRuntime, _require_opencl_binding
from .structs import (
    get_observer_params_struct,
    get_solver_params_struct,
    pack_observer_params,
    pack_solver_params,
)


N_RNGSTATE = 2
@dataclass(slots=True)
class CommonBuffers:
    ensemble_size: int
    problem_shape: ProblemShape
    tspan: object
    solver_params: object
    x0: object
    pars: object
    xf: object
    rng_state: object
    dt: object
    tf: object


@dataclass(slots=True)
class TrajectoryBuffers:
    max_store: int
    t: object
    x: object
    dx: object
    aux: object
    n_stored: object


@dataclass(slots=True)
class FeatureBuffers:
    n_features: int
    observer_data_nbytes: int
    observer_data: object
    observer_params: object
    features: object


class ArrayLayout:
    @staticmethod
    def real_dtype(precision: Precision) -> np.dtype:
        return np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)

    @staticmethod
    def flatten_problem_matrix(
        array: np.ndarray, precision: Precision | None = None
    ) -> np.ndarray:
        dtype = np.float64 if precision is None else ArrayLayout.real_dtype(precision)
        return np.asarray(array, dtype=dtype).flatten(order="F")

    @staticmethod
    def reshape_state(array: np.ndarray | list[float], ensemble_size: int, width: int) -> np.ndarray:
        return np.asarray(array, dtype=np.float64).reshape(
            (ensemble_size, width), order="F"
        )

    @staticmethod
    def reshape_vector(array: np.ndarray | list[float], shape: tuple[int, ...]) -> np.ndarray:
        return np.asarray(array, dtype=np.float64).reshape(shape, order="F")


class BufferManager:
    def __init__(self, runtime: OpenCLRuntime, precision: Precision):
        self._runtime = runtime
        self._precision = precision
        self._real_dtype = ArrayLayout.real_dtype(precision)
        self._opencl_binding = _require_opencl_binding()

    @property
    def real_dtype(self) -> np.dtype:
        return self._real_dtype

    def allocate_common(self, ensemble_size: int, shape: ProblemShape) -> CommonBuffers:
        flags = self._opencl_binding.mem_flags
        real_bytes = self._real_dtype.itemsize
        x0_elements = max(1, ensemble_size * shape.n_var)
        pars_elements = max(1, ensemble_size * shape.n_par)
        rng_elements = max(1, ensemble_size * N_RNGSTATE)

        return CommonBuffers(
            ensemble_size=ensemble_size,
            problem_shape=shape,
            tspan=self._opencl_binding.Buffer(
                self._runtime.context, flags.READ_ONLY, size=2 * real_bytes
            ),
            solver_params=self._opencl_binding.Buffer(
                self._runtime.context,
                flags.READ_ONLY,
                size=get_solver_params_struct(self._runtime, self._precision).dtype.itemsize,
            ),
            x0=self._opencl_binding.Buffer(
                self._runtime.context, flags.READ_WRITE, size=x0_elements * real_bytes
            ),
            pars=self._opencl_binding.Buffer(
                self._runtime.context, flags.READ_ONLY, size=pars_elements * real_bytes
            ),
            xf=self._opencl_binding.Buffer(
                self._runtime.context, flags.READ_WRITE, size=x0_elements * real_bytes
            ),
            rng_state=self._opencl_binding.Buffer(
                self._runtime.context,
                flags.READ_WRITE,
                size=rng_elements * np.dtype(np.uint64).itemsize,
            ),
            dt=self._opencl_binding.Buffer(
                self._runtime.context, flags.READ_WRITE, size=ensemble_size * real_bytes
            ),
            tf=self._opencl_binding.Buffer(
                self._runtime.context, flags.WRITE_ONLY, size=ensemble_size * real_bytes
            ),
        )

    def allocate_trajectory(
        self, ensemble_size: int, shape: ProblemShape, max_store: int
    ) -> TrajectoryBuffers:
        flags = self._opencl_binding.mem_flags
        real_bytes = self._real_dtype.itemsize
        t_elements = max(1, ensemble_size * max_store)
        state_elements = max(1, ensemble_size * shape.n_var * max_store)
        aux_elements = max(1, ensemble_size * shape.n_aux * max_store)

        return TrajectoryBuffers(
            max_store=max_store,
            t=self._opencl_binding.Buffer(
                self._runtime.context, flags.WRITE_ONLY, size=t_elements * real_bytes
            ),
            x=self._opencl_binding.Buffer(
                self._runtime.context, flags.WRITE_ONLY, size=state_elements * real_bytes
            ),
            dx=self._opencl_binding.Buffer(
                self._runtime.context, flags.WRITE_ONLY, size=state_elements * real_bytes
            ),
            aux=self._opencl_binding.Buffer(
                self._runtime.context, flags.WRITE_ONLY, size=aux_elements * real_bytes
            ),
            n_stored=self._opencl_binding.Buffer(
                self._runtime.context,
                flags.WRITE_ONLY,
                size=ensemble_size * np.dtype(np.int32).itemsize,
            ),
        )

    def allocate_feature(
        self,
        ensemble_size: int,
        n_features: int,
        observer_data_nbytes: int,
    ) -> FeatureBuffers:
        flags = self._opencl_binding.mem_flags
        real_bytes = self._real_dtype.itemsize
        feature_elements = max(1, ensemble_size * n_features)
        observer_bytes = max(1, ensemble_size * observer_data_nbytes)
        return FeatureBuffers(
            n_features=n_features,
            observer_data_nbytes=observer_data_nbytes,
            observer_data=self._opencl_binding.Buffer(
                self._runtime.context,
                flags.READ_WRITE,
                size=observer_bytes,
            ),
            observer_params=self._opencl_binding.Buffer(
                self._runtime.context,
                flags.READ_ONLY,
                size=get_observer_params_struct(
                    self._runtime, self._precision
                ).dtype.itemsize,
            ),
            features=self._opencl_binding.Buffer(
                self._runtime.context,
                flags.WRITE_ONLY,
                size=feature_elements * real_bytes,
            ),
        )

    def upload_tspan(
        self, buffers: CommonBuffers, tspan: tuple[float, float]
    ) -> np.ndarray:
        host = np.asarray(tspan, dtype=self._real_dtype)
        self._enqueue_copy(buffers.tspan, host)
        return host

    def upload_solver_params(
        self, buffers: CommonBuffers, solver_params: SolverParams
    ) -> np.ndarray:
        host = pack_solver_params(self._runtime, solver_params, self._precision)
        self._enqueue_copy(buffers.solver_params, host)
        return host

    def upload_observer_params(
        self, buffers: FeatureBuffers, observer_params: ObserverParams
    ) -> np.ndarray:
        host = pack_observer_params(self._runtime, observer_params, self._precision)
        self._enqueue_copy(buffers.observer_params, host)
        return host

    def clear_observer_data(self, buffers: FeatureBuffers, ensemble_size: int) -> None:
        host = np.zeros(
            max(1, ensemble_size * buffers.observer_data_nbytes),
            dtype=np.uint8,
        )
        self._enqueue_copy(buffers.observer_data, host)

    def upload_problem_data(
        self,
        buffers: CommonBuffers,
        initial_state: np.ndarray,
        parameters: np.ndarray,
    ) -> tuple[np.ndarray, np.ndarray]:
        initial_state_host = self.upload_x0(buffers, initial_state)
        parameters_host = self.upload_pars(buffers, parameters)
        return initial_state_host, parameters_host

    def upload_x0(self, buffers: CommonBuffers, initial_state: np.ndarray) -> np.ndarray:
        host = ArrayLayout.flatten_problem_matrix(initial_state, precision=self._precision)
        self._enqueue_copy(buffers.x0, host)
        return host

    def upload_pars(self, buffers: CommonBuffers, parameters: np.ndarray) -> np.ndarray:
        host = ArrayLayout.flatten_problem_matrix(parameters, precision=self._precision)
        if host.size > 0:
            self._enqueue_copy(buffers.pars, host)
        return host

    def upload_dt(self, buffers: CommonBuffers, dt_values: np.ndarray) -> np.ndarray:
        host = np.asarray(dt_values, dtype=self._real_dtype)
        self._enqueue_copy(buffers.dt, host)
        return host

    def upload_rng_state(
        self, buffers: CommonBuffers, rng_state: np.ndarray
    ) -> np.ndarray:
        host = np.asarray(rng_state, dtype=np.uint64).flatten(order="F")
        self._enqueue_copy(buffers.rng_state, host)
        return host

    def download_x0(self, buffers: CommonBuffers) -> np.ndarray:
        return self._download_state_buffer(buffers.x0, buffers.ensemble_size, buffers.problem_shape.n_var)

    def download_pars(self, buffers: CommonBuffers) -> np.ndarray:
        return self._download_state_buffer(buffers.pars, buffers.ensemble_size, buffers.problem_shape.n_par)

    def download_xf(self, buffers: CommonBuffers) -> np.ndarray:
        return self._download_state_buffer(buffers.xf, buffers.ensemble_size, buffers.problem_shape.n_var)

    def download_dt(self, buffers: CommonBuffers) -> np.ndarray:
        return self._download_vector_buffer(buffers.dt, (buffers.ensemble_size,))

    def download_tf(self, buffers: CommonBuffers) -> np.ndarray:
        return self._download_vector_buffer(buffers.tf, (buffers.ensemble_size,))

    def download_t(
        self, buffers: CommonBuffers, trajectory_buffers: TrajectoryBuffers
    ) -> np.ndarray:
        return self._download_vector_buffer(
            trajectory_buffers.t,
            (buffers.ensemble_size * trajectory_buffers.max_store,),
        )

    def download_x(
        self, buffers: CommonBuffers, trajectory_buffers: TrajectoryBuffers
    ) -> np.ndarray:
        return self._download_vector_buffer(
            trajectory_buffers.x,
            (
                buffers.ensemble_size
                * buffers.problem_shape.n_var
                * trajectory_buffers.max_store,
            ),
        )

    def download_dx(
        self, buffers: CommonBuffers, trajectory_buffers: TrajectoryBuffers
    ) -> np.ndarray:
        return self._download_vector_buffer(
            trajectory_buffers.dx,
            (
                buffers.ensemble_size
                * buffers.problem_shape.n_var
                * trajectory_buffers.max_store,
            ),
        )

    def download_aux(
        self, buffers: CommonBuffers, trajectory_buffers: TrajectoryBuffers
    ) -> np.ndarray:
        if buffers.problem_shape.n_aux == 0:
            return self._download_vector_buffer(trajectory_buffers.aux, (1,))
        return self._download_vector_buffer(
            trajectory_buffers.aux,
            (
                buffers.ensemble_size
                * buffers.problem_shape.n_aux
                * trajectory_buffers.max_store,
            ),
        )

    def download_n_stored(
        self, buffers: CommonBuffers, trajectory_buffers: TrajectoryBuffers
    ) -> np.ndarray:
        host = np.empty(buffers.ensemble_size, dtype=np.int32)
        self._opencl_binding.enqueue_copy(
            self._runtime.queue, host, trajectory_buffers.n_stored, is_blocking=True
        )
        return host

    def download_features(
        self, buffers: CommonBuffers, feature_buffers: FeatureBuffers
    ) -> np.ndarray:
        return self._download_vector_buffer(
            feature_buffers.features,
            (buffers.ensemble_size * feature_buffers.n_features,),
        )

    def download_rng_state(self, buffers: CommonBuffers) -> np.ndarray:
        host = np.empty(buffers.ensemble_size * N_RNGSTATE, dtype=np.uint64)
        self._opencl_binding.enqueue_copy(
            self._runtime.queue, host, buffers.rng_state, is_blocking=True
        )
        return host.reshape((buffers.ensemble_size, N_RNGSTATE), order="F")

    def _download_state_buffer(
        self, buffer: object, ensemble_size: int, width: int
    ) -> np.ndarray:
        if width == 0:
            return np.empty((ensemble_size, 0), dtype=np.float64)
        host = np.empty(ensemble_size * width, dtype=self._real_dtype)
        self._opencl_binding.enqueue_copy(
            self._runtime.queue, host, buffer, is_blocking=True
        )
        return ArrayLayout.reshape_state(host, ensemble_size, width)

    def _download_vector_buffer(
        self, buffer: object, shape: tuple[int, ...]
    ) -> np.ndarray:
        host = np.empty(int(np.prod(shape)), dtype=self._real_dtype)
        self._opencl_binding.enqueue_copy(
            self._runtime.queue, host, buffer, is_blocking=True
        )
        return ArrayLayout.reshape_vector(host, shape)

    def _enqueue_copy(self, buffer: object, host: np.ndarray) -> None:
        if host.size == 0:
            return
        self._opencl_binding.enqueue_copy(
            self._runtime.queue, buffer, host, is_blocking=True
        )
