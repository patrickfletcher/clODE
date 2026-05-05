from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

if TYPE_CHECKING:
    from clode.cpp.clode_cpp_wrapper import (
        ObserverParams as CppObserverParams,
        ProblemInfo as CppProblemInfo,
        SolverParams as CppSolverParams,
    )


@dataclass(slots=True)
class ProblemInfo:
    src_file: str = ""
    vars: list[str] = field(default_factory=list)
    pars: list[str] = field(default_factory=list)
    aux: list[str] = field(default_factory=list)
    num_noise: int = 1

    def __post_init__(self) -> None:
        self.vars = list(self.vars)
        self.pars = list(self.pars)
        self.aux = list(self.aux)
        self.num_noise = int(self.num_noise)
        if self.num_noise < 0:
            raise ValueError("num_noise must be non-negative")

    @property
    def num_var(self) -> int:
        return len(self.vars)

    @property
    def num_par(self) -> int:
        return len(self.pars)

    @property
    def num_aux(self) -> int:
        return len(self.aux)


@dataclass(slots=True)
class SolverParams:
    dt: float = 0.1
    dtmax: float = 0.5
    abstol: float = 1e-6
    reltol: float = 1e-3
    max_steps: int = 1_000_000
    max_store: int = 1_000_000
    nout: int = 1

    def __post_init__(self) -> None:
        self.dt = float(self.dt)
        self.dtmax = float(self.dtmax)
        self.abstol = float(self.abstol)
        self.reltol = float(self.reltol)
        self.max_steps = int(self.max_steps)
        self.max_store = int(self.max_store)
        self.nout = int(self.nout)


@dataclass(slots=True)
class ObserverParams:
    e_var_ix: int = 0
    f_var_ix: int = 0
    max_event_count: int = 100
    max_event_timestamps: int = 0
    min_amp: float = 0.0
    min_imi: float = 0.0
    nhood_radius: float = 0.05
    x_up_threshold: float = 0.2
    x_down_threshold: float = 0.2
    dx_up_threshold: float = 0.0
    dx_down_threshold: float = 0.0
    eps_dx: float = 0.0

    def __post_init__(self) -> None:
        self.e_var_ix = int(self.e_var_ix)
        self.f_var_ix = int(self.f_var_ix)
        self.max_event_count = int(self.max_event_count)
        self.max_event_timestamps = int(self.max_event_timestamps)
        self.min_amp = float(self.min_amp)
        self.min_imi = float(self.min_imi)
        self.nhood_radius = float(self.nhood_radius)
        self.x_up_threshold = float(self.x_up_threshold)
        self.x_down_threshold = float(self.x_down_threshold)
        self.dx_up_threshold = float(self.dx_up_threshold)
        self.dx_down_threshold = float(self.dx_down_threshold)
        self.eps_dx = float(self.eps_dx)


def problem_info_from_cpp(problem_info: Any) -> ProblemInfo:
    return ProblemInfo(
        src_file=str(problem_info.src_file),
        vars=list(problem_info.vars),
        pars=list(problem_info.pars),
        aux=list(problem_info.aux),
        num_noise=int(problem_info.num_noise),
    )


def solver_params_from_cpp(solver_params: Any) -> SolverParams:
    return SolverParams(
        dt=solver_params.dt,
        dtmax=solver_params.dtmax,
        abstol=solver_params.abstol,
        reltol=solver_params.reltol,
        max_steps=solver_params.max_steps,
        max_store=solver_params.max_store,
        nout=solver_params.nout,
    )


def observer_params_from_cpp(observer_params: Any) -> ObserverParams:
    return ObserverParams(
        e_var_ix=observer_params.e_var_ix,
        f_var_ix=observer_params.f_var_ix,
        max_event_count=observer_params.max_event_count,
        max_event_timestamps=observer_params.max_event_timestamps,
        min_amp=observer_params.min_amp,
        min_imi=observer_params.min_imi,
        nhood_radius=observer_params.nhood_radius,
        x_up_threshold=observer_params.x_up_threshold,
        x_down_threshold=observer_params.x_down_threshold,
        dx_up_threshold=observer_params.dx_up_threshold,
        dx_down_threshold=observer_params.dx_down_threshold,
        eps_dx=observer_params.eps_dx,
    )


def problem_info_to_cpp(problem_info: ProblemInfo) -> CppProblemInfo:
    from clode.cpp.clode_cpp_wrapper import ProblemInfo as CppProblemInfo

    return CppProblemInfo(
        problem_info.src_file,
        list(problem_info.vars),
        list(problem_info.pars),
        list(problem_info.aux),
        problem_info.num_noise,
    )


def solver_params_to_cpp(solver_params: SolverParams) -> CppSolverParams:
    from clode.cpp.clode_cpp_wrapper import SolverParams as CppSolverParams

    return CppSolverParams(
        solver_params.dt,
        solver_params.dtmax,
        solver_params.abstol,
        solver_params.reltol,
        solver_params.max_steps,
        solver_params.max_store,
        solver_params.nout,
    )


def observer_params_to_cpp(observer_params: ObserverParams) -> CppObserverParams:
    from clode.cpp.clode_cpp_wrapper import ObserverParams as CppObserverParams

    return CppObserverParams(
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


__all__ = [
    "ObserverParams",
    "ProblemInfo",
    "SolverParams",
    "observer_params_from_cpp",
    "observer_params_to_cpp",
    "problem_info_from_cpp",
    "problem_info_to_cpp",
    "solver_params_from_cpp",
    "solver_params_to_cpp",
]