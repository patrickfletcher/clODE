from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class _TrajectoryOutputSettings:
    max_store: int = 1_000_000
    nout: int = 1

    def __post_init__(self) -> None:
        object.__setattr__(self, "max_store", int(self.max_store))
        object.__setattr__(self, "nout", int(self.nout))


@dataclass(frozen=True, slots=True)
class _IntegrationSettings:
    dt: float = 0.1
    dtmax: float = 0.5
    abstol: float = 1e-6
    reltol: float = 1e-3
    max_steps: int = 1_000_000

    def __post_init__(self) -> None:
        object.__setattr__(self, "dt", float(self.dt))
        object.__setattr__(self, "dtmax", float(self.dtmax))
        object.__setattr__(self, "abstol", float(self.abstol))
        object.__setattr__(self, "reltol", float(self.reltol))
        object.__setattr__(self, "max_steps", int(self.max_steps))

    def to_solver_params(
        self, output_settings: _TrajectoryOutputSettings
    ) -> SolverParams:
        return SolverParams(
            dt=self.dt,
            dtmax=self.dtmax,
            abstol=self.abstol,
            reltol=self.reltol,
            max_steps=self.max_steps,
            max_store=output_settings.max_store,
            nout=output_settings.nout,
        )


@dataclass(slots=True)
class SolverParams:
    """Compatibility bundle for integration and trajectory-output settings.

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

    def copy(self) -> SolverParams:
        return SolverParams(
            dt=self.dt,
            dtmax=self.dtmax,
            abstol=self.abstol,
            reltol=self.reltol,
            max_steps=self.max_steps,
            max_store=self.max_store,
            nout=self.nout,
        )

    @property
    def integration_settings(self) -> _IntegrationSettings:
        return _IntegrationSettings(
            dt=self.dt,
            dtmax=self.dtmax,
            abstol=self.abstol,
            reltol=self.reltol,
            max_steps=self.max_steps,
        )

    @property
    def trajectory_output_settings(self) -> _TrajectoryOutputSettings:
        return _TrajectoryOutputSettings(
            max_store=self.max_store,
            nout=self.nout,
        )


__all__ = ["SolverParams"]