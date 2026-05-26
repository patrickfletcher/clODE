from __future__ import annotations

from dataclasses import dataclass


_MAX_UINT64 = (1 << 64) - 1


def _coerce_max_steps(value: int) -> int:
    max_steps = int(value)
    if max_steps < 0:
        raise ValueError("max_steps must be non-negative")
    if max_steps > _MAX_UINT64:
        raise OverflowError("max_steps exceeds uint64 range")
    return max_steps


@dataclass(frozen=True, slots=True)
class TrajectoryOutputSettings:
    """Internal value object for retained trajectory-sample policy."""

    max_store: int = 1_000_000
    nout: int = 1

    def __post_init__(self) -> None:
        object.__setattr__(self, "max_store", int(self.max_store))
        object.__setattr__(self, "nout", int(self.nout))


@dataclass(frozen=True, slots=True)
class IntegrationSettings:
    """Internal value object for integration policy independent of output policy."""

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
        object.__setattr__(self, "max_steps", _coerce_max_steps(self.max_steps))

    def to_solver_params(
        self, output_settings: TrajectoryOutputSettings
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


_DEFAULT_TRAJECTORY_OUTPUT_SETTINGS = TrajectoryOutputSettings()
_DEFAULT_INTEGRATION_SETTINGS = IntegrationSettings()


@dataclass(slots=True)
class SolverParams:
    """Public compatibility bundle for integration and trajectory-output settings.

    This object is retained for the current user-facing API. Internal code should
    prefer `IntegrationSettings` and `TrajectoryOutputSettings` when the narrower
    owner model is sufficient.

    Attributes:
        dt: Initial or fixed time step.
        dtmax: Maximum time step for adaptive steppers.
        abstol: Absolute tolerance for adaptive steppers.
        reltol: Relative tolerance for adaptive steppers.
        max_steps: Maximum number of integration steps per solve.
        max_store: Maximum number of stored trajectory samples.
        nout: Storage/output stride used by the trajectory path.
    """

    dt: float = _DEFAULT_INTEGRATION_SETTINGS.dt
    dtmax: float = _DEFAULT_INTEGRATION_SETTINGS.dtmax
    abstol: float = _DEFAULT_INTEGRATION_SETTINGS.abstol
    reltol: float = _DEFAULT_INTEGRATION_SETTINGS.reltol
    max_steps: int = _DEFAULT_INTEGRATION_SETTINGS.max_steps
    max_store: int = _DEFAULT_TRAJECTORY_OUTPUT_SETTINGS.max_store
    nout: int = _DEFAULT_TRAJECTORY_OUTPUT_SETTINGS.nout

    def __post_init__(self) -> None:
        self.dt = float(self.dt)
        self.dtmax = float(self.dtmax)
        self.abstol = float(self.abstol)
        self.reltol = float(self.reltol)
        self.max_steps = _coerce_max_steps(self.max_steps)
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
    def integration_settings(self) -> IntegrationSettings:
        return IntegrationSettings(
            dt=self.dt,
            dtmax=self.dtmax,
            abstol=self.abstol,
            reltol=self.reltol,
            max_steps=self.max_steps,
        )

    @property
    def trajectory_output_settings(self) -> TrajectoryOutputSettings:
        return TrajectoryOutputSettings(
            max_store=self.max_store,
            nout=self.nout,
        )


def _resolve_solver_params(
    *,
    solver_parameters: SolverParams | None = None,
    dt: float = _DEFAULT_INTEGRATION_SETTINGS.dt,
    dtmax: float = _DEFAULT_INTEGRATION_SETTINGS.dtmax,
    abstol: float = _DEFAULT_INTEGRATION_SETTINGS.abstol,
    reltol: float = _DEFAULT_INTEGRATION_SETTINGS.reltol,
    max_steps: int = _DEFAULT_INTEGRATION_SETTINGS.max_steps,
    max_store: int = _DEFAULT_TRAJECTORY_OUTPUT_SETTINGS.max_store,
    nout: int = _DEFAULT_TRAJECTORY_OUTPUT_SETTINGS.nout,
) -> SolverParams:
    if solver_parameters is not None:
        return solver_parameters.copy()

    integration_settings = IntegrationSettings(
        dt=dt,
        dtmax=dtmax,
        abstol=abstol,
        reltol=reltol,
        max_steps=max_steps,
    )
    output_settings = TrajectoryOutputSettings(
        max_store=max_store,
        nout=nout,
    )
    return integration_settings.to_solver_params(output_settings)


__all__ = ["SolverParams"]