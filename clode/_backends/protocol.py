from __future__ import annotations

from typing import Protocol, Sequence

from ..observers.types import ObserverParams
from ..simulation.params import SolverParams


class SimulatorBackend(Protocol):
    def build_cl(self) -> None:
        ...

    def get_available_steppers(self) -> list[str]:
        ...

    def get_dt(self) -> list[float]:
        ...

    def get_program_string(self) -> str:
        ...

    def get_solver_params(self) -> SolverParams:
        ...

    def get_tf(self) -> list[float]:
        ...

    def get_tspan(self) -> list[float]:
        ...

    def get_x0(self) -> list[float]:
        ...

    def get_xf(self) -> list[float]:
        ...

    def print_status(self) -> None:
        ...

    def seed_rng(self, seed: int | None = None) -> None:
        ...

    def set_pars(self, parameters: Sequence[float]) -> None:
        ...

    def set_problem_data(
        self, initial_state: Sequence[float], parameters: Sequence[float]
    ) -> None:
        ...

    def set_solver_params(self, solver_params: SolverParams) -> None:
        ...

    def set_tspan(self, tspan: Sequence[float]) -> None:
        ...

    def set_x0(self, initial_state: Sequence[float]) -> None:
        ...

    def shift_tspan(self) -> None:
        ...

    def shift_x0(self) -> None:
        ...

    def transient(self) -> None:
        ...


class TrajectoryBackend(SimulatorBackend, Protocol):
    def get_aux(self) -> list[float]:
        ...

    def get_dx(self) -> list[float]:
        ...

    def get_n_stored(self) -> list[int]:
        ...

    def get_t(self) -> list[float]:
        ...

    def get_x(self) -> list[float]:
        ...

    def trajectory(self) -> None:
        ...


class FeatureBackend(SimulatorBackend, Protocol):
    def features(self, reinitialize_observer: bool | None = None) -> None:
        ...

    def get_f(self) -> list[float]:
        ...

    def get_feature_names(self) -> list[str]:
        ...

    def get_n_features(self) -> int:
        ...

    def get_observer_params(self) -> ObserverParams:
        ...

    def initialize_observer(self) -> None:
        ...

    def is_observer_initialized(self) -> bool:
        ...

    def set_observer(self, observer: str) -> None:
        ...

    def set_observer_params(self, observer_params: ObserverParams) -> None:
        ...
