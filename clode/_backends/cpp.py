from __future__ import annotations

from typing import Generic, Sequence, TypeVar

from clode.cpp.clode_cpp_wrapper import (
    FeatureSimulatorBase,
    ObserverParams,
    ProblemInfo,
    SimulatorBase,
    SolverParams,
    TrajectorySimulatorBase,
)

from ..runtime import OpenCLResource

PybindBackend = TypeVar(
    "PybindBackend", SimulatorBase, TrajectorySimulatorBase, FeatureSimulatorBase
)


class _CppBackendBase(Generic[PybindBackend]):
    def __init__(self, integrator: PybindBackend) -> None:
        self._integrator = integrator

    def build_cl(self) -> None:
        self._integrator.build_cl()

    def get_available_steppers(self) -> list[str]:
        return self._integrator.get_available_steppers()

    def get_dt(self) -> list[float]:
        return self._integrator.get_dt()

    def get_program_string(self) -> str:
        return self._integrator.get_program_string()

    def get_solver_params(self) -> SolverParams:
        return self._integrator.get_solver_params()

    def get_tf(self) -> list[float]:
        return self._integrator.get_tf()

    def get_tspan(self) -> list[float]:
        return self._integrator.get_tspan()

    def get_x0(self) -> list[float]:
        return self._integrator.get_x0()

    def get_xf(self) -> list[float]:
        return self._integrator.get_xf()

    def print_status(self) -> None:
        self._integrator.print_status()

    def seed_rng(self, seed: int | None = None) -> None:
        if seed is None:
            self._integrator.seed_rng()
        else:
            self._integrator.seed_rng(seed)

    def set_pars(self, parameters: Sequence[float]) -> None:
        self._integrator.set_pars(list(parameters))

    def set_problem_data(
        self, initial_state: Sequence[float], parameters: Sequence[float]
    ) -> None:
        self._integrator.set_problem_data(list(initial_state), list(parameters))

    def set_solver_params(self, solver_params: SolverParams) -> None:
        self._integrator.set_solver_params(solver_params)

    def set_tspan(self, tspan: Sequence[float]) -> None:
        self._integrator.set_tspan(list(tspan))

    def set_x0(self, initial_state: Sequence[float]) -> None:
        self._integrator.set_x0(list(initial_state))

    def shift_tspan(self) -> None:
        self._integrator.shift_tspan()

    def shift_x0(self) -> None:
        self._integrator.shift_x0()

    def transient(self) -> None:
        self._integrator.transient()


class CppSimulatorBackend(_CppBackendBase[SimulatorBase]):
    def __init__(
        self,
        problem_info: ProblemInfo,
        stepper: str,
        single_precision: bool,
        runtime: OpenCLResource,
        clode_root: str,
    ) -> None:
        super().__init__(
            SimulatorBase(problem_info, stepper, single_precision, runtime, clode_root)
        )


class CppTrajectoryBackend(_CppBackendBase[TrajectorySimulatorBase]):
    def __init__(
        self,
        problem_info: ProblemInfo,
        stepper: str,
        single_precision: bool,
        runtime: OpenCLResource,
        clode_root: str,
    ) -> None:
        super().__init__(
            TrajectorySimulatorBase(
                problem_info, stepper, single_precision, runtime, clode_root
            )
        )

    def get_aux(self) -> list[float]:
        return self._integrator.get_aux()

    def get_dx(self) -> list[float]:
        return self._integrator.get_dx()

    def get_n_stored(self) -> list[int]:
        return self._integrator.get_n_stored()

    def get_t(self) -> list[float]:
        return self._integrator.get_t()

    def get_x(self) -> list[float]:
        return self._integrator.get_x()

    def trajectory(self) -> None:
        self._integrator.trajectory()


class CppFeatureBackend(_CppBackendBase[FeatureSimulatorBase]):
    def __init__(
        self,
        problem_info: ProblemInfo,
        stepper: str,
        observer: str,
        observer_params: ObserverParams,
        single_precision: bool,
        runtime: OpenCLResource,
        clode_root: str,
    ) -> None:
        super().__init__(
            FeatureSimulatorBase(
                problem_info,
                stepper,
                observer,
                observer_params,
                single_precision,
                runtime,
                clode_root,
            )
        )

    def features(self, reinitialize_observer: bool | None = None) -> None:
        if reinitialize_observer is None:
            self._integrator.features()
        else:
            self._integrator.features(reinitialize_observer)

    def get_f(self) -> list[float]:
        return self._integrator.get_f()

    def get_feature_names(self) -> list[str]:
        return self._integrator.get_feature_names()

    def get_n_features(self) -> int:
        return self._integrator.get_n_features()

    def get_observer_params(self) -> ObserverParams:
        return self._integrator.get_observer_params()

    def initialize_observer(self) -> None:
        self._integrator.initialize_observer()

    def is_observer_initialized(self) -> bool:
        return self._integrator.is_observer_initialized()

    def set_observer(self, observer: str) -> None:
        self._integrator.set_observer(observer)

    def set_observer_params(self, observer_params: ObserverParams) -> None:
        self._integrator.set_observer_params(observer_params)
