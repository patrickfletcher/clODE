from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
from typing import Any

import numpy as np

FloatArray = np.ndarray[Any, np.dtype[np.float64]]
IntArray = np.ndarray[Any, np.dtype[np.int32]]
UIntArray = np.ndarray[Any, np.dtype[np.uint64]]


class SolverStatus(IntEnum):
    """Per-instance solve completion or failure status.

    Keep these values in sync with the kernel-side status codes written by the
    transient, feature, and trajectory entrypoints.
    """

    COMPLETED = 0
    MAX_STEPS_REACHED = 1
    TERMINAL_EVENT_REACHED = 2
    OUTPUT_CAPACITY_REACHED = 3
    NO_PROGRESS = 4
    STEPPER_FAILED = -1


@dataclass(slots=True)
class SolverState:
    t_span: tuple[float, float] = (0.0, 0.0)
    current_time: FloatArray | None = None
    current_dt: FloatArray | None = None
    status: IntArray | None = None
    step_count: UIntArray | None = None
    last_accepted_dt: FloatArray | None = None
    final_time: FloatArray | None = None
    problem_data_needs_pull: bool = False

    def set_requested_window(self, t_span: tuple[float, float]) -> None:
        self.t_span = (float(t_span[0]), float(t_span[1]))

    def reset_problem_time(self, ensemble_shape: tuple[int, ...]) -> None:
        self.current_time = np.full(ensemble_shape, self.t_span[0], dtype=np.float64)
        self.problem_data_needs_pull = False

    def continue_problem_time(self) -> None:
        self.current_time = (
            None
            if self.final_time is None
            else np.array(self.final_time, dtype=np.float64, copy=True)
        )
        self.problem_data_needs_pull = True

    def attained_final_time_window(
        self,
        *,
        atol: float = 1e-12,
        rtol: float = 0.0,
    ) -> tuple[float, float]:
        if self.final_time is None:
            raise ValueError("Must run a simulation before getting final time")

        final_times = np.asarray(self.final_time, dtype=np.float64).reshape(-1)
        if final_times.size == 0:
            raise ValueError("Must run a simulation before getting final time")
        if not np.allclose(final_times, final_times[0], atol=atol, rtol=rtol):
            raise ValueError(
                "Cannot advance to one shared continuation window when ensemble members reached different final times"
            )

        next_start = float(final_times[0])
        duration = float(self.t_span[1] - self.t_span[0])
        return (next_start, next_start + duration)

    def mark_problem_data_synced(self) -> None:
        self.problem_data_needs_pull = False

    def invalidate_results(self) -> None:
        self.current_dt = None
        self.status = None
        self.step_count = None
        self.last_accepted_dt = None
        self.final_time = None


@dataclass(slots=True)
class TransientCache:
    final_state: FloatArray | None = None

    def invalidate(self) -> None:
        self.final_state = None


@dataclass(slots=True)
class FeatureCache:
    feature_array: FloatArray | None = None
    num_features: int | None = None
    has_result: bool = False

    def invalidate(self) -> None:
        self.feature_array = None
        self.num_features = None
        self.has_result = False

    def mark_result_pending(self) -> None:
        self.feature_array = None
        self.num_features = None
        self.has_result = True


@dataclass(slots=True)
class TrajectoryCache:
    n_stored: IntArray | None = None
    t: FloatArray | None = None
    x: FloatArray | None = None
    dx: FloatArray | None = None
    aux: FloatArray | None = None
    has_result: bool = False

    def invalidate(self) -> None:
        self.n_stored = None
        self.t = None
        self.x = None
        self.dx = None
        self.aux = None
        self.has_result = False

    def mark_result_pending(self) -> None:
        self.n_stored = None
        self.t = None
        self.x = None
        self.dx = None
        self.aux = None
        self.has_result = True


__all__ = [
    "FeatureCache",
    "SolverStatus",
    "SolverState",
    "TrajectoryCache",
    "TransientCache",
]