from __future__ import annotations

from dataclasses import dataclass, field


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
__all__ = [
    "ObserverParams",
    "ProblemInfo",
    "SolverParams",
]