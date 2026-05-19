from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

FloatArray = np.ndarray[Any, np.dtype[np.float64]]
IntArray = np.ndarray[Any, np.dtype[np.int32]]


@dataclass(slots=True)
class SolverState:
    t_span: tuple[float, float] = (0.0, 0.0)
    current_time: FloatArray | None = None
    current_dt: FloatArray | None = None
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

    def mark_problem_data_synced(self) -> None:
        self.problem_data_needs_pull = False

    def invalidate_results(self) -> None:
        self.current_dt = None
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
    "SolverState",
    "TrajectoryCache",
    "TransientCache",
]