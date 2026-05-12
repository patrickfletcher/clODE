from __future__ import annotations

from dataclasses import dataclass


@dataclass(slots=True)
class SolverParams:
    """Solver configuration shared by the simulator classes.

    Attributes:
        dt: Initial or fixed time step.
        dtmax: Maximum time step for adaptive steppers.
        abstol: Absolute tolerance for adaptive steppers.
        reltol: Relative tolerance for adaptive steppers.
        max_steps: Maximum number of integration steps per solve.
        max_store: Maximum number of stored trajectory samples.
        nout: Storage/output stride used by the trajectory path.
    """

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


__all__ = ["SolverParams"]