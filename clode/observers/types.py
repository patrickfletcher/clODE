from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum


_DEFAULT_EVENT_VAR_INDEX = 0
_DEFAULT_FEATURE_VAR_INDEX = 0
_DEFAULT_MAX_EVENT_COUNT = 100
_DEFAULT_MAX_EVENT_TIMESTAMPS = 0
_DEFAULT_MIN_AMP = 0.0
_DEFAULT_MIN_IMI = 0.0
_DEFAULT_NHOOD_RADIUS = 0.05
_DEFAULT_X_UP_THRESHOLD = 0.3
_DEFAULT_X_DOWN_THRESHOLD = 0.2
_DEFAULT_DX_UP_THRESHOLD = 0.0
_DEFAULT_DX_DOWN_THRESHOLD = 0.0
_DEFAULT_EPS_DX = 0.0


@dataclass(frozen=True, slots=True)
class ObserverRuntimeSettings:
    """Internal value object for observer runtime thresholds and selector indices."""

    e_var_ix: int = _DEFAULT_EVENT_VAR_INDEX
    f_var_ix: int = _DEFAULT_FEATURE_VAR_INDEX
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    min_amp: float = _DEFAULT_MIN_AMP
    min_imi: float = _DEFAULT_MIN_IMI
    nhood_radius: float = _DEFAULT_NHOOD_RADIUS
    x_up_threshold: float = _DEFAULT_X_UP_THRESHOLD
    x_down_threshold: float = _DEFAULT_X_DOWN_THRESHOLD
    dx_up_threshold: float = _DEFAULT_DX_UP_THRESHOLD
    dx_down_threshold: float = _DEFAULT_DX_DOWN_THRESHOLD
    eps_dx: float = _DEFAULT_EPS_DX

    def __post_init__(self) -> None:
        object.__setattr__(self, "e_var_ix", int(self.e_var_ix))
        object.__setattr__(self, "f_var_ix", int(self.f_var_ix))
        object.__setattr__(self, "max_event_count", int(self.max_event_count))
        object.__setattr__(self, "min_amp", float(self.min_amp))
        object.__setattr__(self, "min_imi", float(self.min_imi))
        object.__setattr__(self, "nhood_radius", float(self.nhood_radius))
        object.__setattr__(self, "x_up_threshold", float(self.x_up_threshold))
        object.__setattr__(self, "x_down_threshold", float(self.x_down_threshold))
        object.__setattr__(self, "dx_up_threshold", float(self.dx_up_threshold))
        object.__setattr__(self, "dx_down_threshold", float(self.dx_down_threshold))
        object.__setattr__(self, "eps_dx", float(self.eps_dx))

    def to_observer_params(
        self, event_output_settings: EventOutputSettings
    ) -> ObserverParams:
        return ObserverParams(
            e_var_ix=self.e_var_ix,
            f_var_ix=self.f_var_ix,
            max_event_count=self.max_event_count,
            max_event_timestamps=event_output_settings.max_event_timestamps,
            min_amp=self.min_amp,
            min_imi=self.min_imi,
            nhood_radius=self.nhood_radius,
            x_up_threshold=self.x_up_threshold,
            x_down_threshold=self.x_down_threshold,
            dx_up_threshold=self.dx_up_threshold,
            dx_down_threshold=self.dx_down_threshold,
            eps_dx=self.eps_dx,
        )


@dataclass(frozen=True, slots=True)
class EventOutputSettings:
    """Internal value object for retained observer event-output policy."""

    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS

    def __post_init__(self) -> None:
        object.__setattr__(self, "max_event_timestamps", int(self.max_event_timestamps))


class Observer(Enum):
    """Built-in observer modes available to `FeatureSimulator`."""

    basic = "basic"
    basic_all_variables = "basicall"
    local_max = "localmax"
    neighbourhood_1 = "nhood1"
    neighbourhood_2 = "nhood2"
    threshold_2 = "thresh2"


@dataclass(slots=True)
class ObserverParams:
    """Public compatibility bundle for built-in observer configuration.

    This object is retained for the current user-facing API. Internal code should
    prefer `ObserverRuntimeSettings` and `EventOutputSettings` when the narrower
    owner model is sufficient.

    Attributes:
        e_var_ix: Index of the variable used for event detection.
        f_var_ix: Index of the variable used for feature readout.
        max_event_count: Maximum number of events to accumulate.
        max_event_timestamps: Maximum number of event timestamps to retain.
        min_amp: Minimum amplitude threshold for event acceptance.
        min_imi: Minimum inter-event interval.
        nhood_radius: Neighborhood radius used by neighborhood observers.
        x_up_threshold: Rising threshold for threshold-style observers.
        x_down_threshold: Falling threshold for threshold-style observers.
        dx_up_threshold: Rising derivative threshold for threshold-style observers.
        dx_down_threshold: Falling derivative threshold for threshold-style observers.
        eps_dx: Small derivative tolerance used around threshold crossings.
    """

    e_var_ix: int = _DEFAULT_EVENT_VAR_INDEX
    f_var_ix: int = _DEFAULT_FEATURE_VAR_INDEX
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS
    min_amp: float = _DEFAULT_MIN_AMP
    min_imi: float = _DEFAULT_MIN_IMI
    nhood_radius: float = _DEFAULT_NHOOD_RADIUS
    x_up_threshold: float = _DEFAULT_X_UP_THRESHOLD
    x_down_threshold: float = _DEFAULT_X_DOWN_THRESHOLD
    dx_up_threshold: float = _DEFAULT_DX_UP_THRESHOLD
    dx_down_threshold: float = _DEFAULT_DX_DOWN_THRESHOLD
    eps_dx: float = _DEFAULT_EPS_DX

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

    @property
    def runtime_settings(self) -> ObserverRuntimeSettings:
        return ObserverRuntimeSettings(
            e_var_ix=self.e_var_ix,
            f_var_ix=self.f_var_ix,
            max_event_count=self.max_event_count,
            min_amp=self.min_amp,
            min_imi=self.min_imi,
            nhood_radius=self.nhood_radius,
            x_up_threshold=self.x_up_threshold,
            x_down_threshold=self.x_down_threshold,
            dx_up_threshold=self.dx_up_threshold,
            dx_down_threshold=self.dx_down_threshold,
            eps_dx=self.eps_dx,
        )

    @property
    def event_output_settings(self) -> EventOutputSettings:
        return EventOutputSettings(
            max_event_timestamps=self.max_event_timestamps,
        )


def _copy_observer_params(observer_params: ObserverParams) -> ObserverParams:
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


def _resolve_observer_params(
    variable_names: Sequence[str],
    *,
    observer_params: ObserverParams | None = None,
    base_params: ObserverParams | None = None,
    event_var: str | None = None,
    feature_var: str | None = None,
    max_event_count: int | None = None,
    max_event_timestamps: int | None = None,
    min_amp: float | None = None,
    min_imi: float | None = None,
    nhood_radius: float | None = None,
    x_up_threshold: float | None = None,
    x_down_threshold: float | None = None,
    dx_up_threshold: float | None = None,
    dx_down_threshold: float | None = None,
    eps_dx: float | None = None,
) -> ObserverParams:
    if observer_params is not None:
        return _copy_observer_params(observer_params)

    current = ObserverParams() if base_params is None else _copy_observer_params(base_params)
    return ObserverParams(
        e_var_ix=_resolve_variable_index(
            variable_names,
            event_var,
            default_index=current.e_var_ix,
            parameter_name="event_var",
        ),
        f_var_ix=_resolve_variable_index(
            variable_names,
            feature_var,
            default_index=current.f_var_ix,
            parameter_name="feature_var",
        ),
        max_event_count=current.max_event_count if max_event_count is None else max_event_count,
        max_event_timestamps=(
            current.max_event_timestamps
            if max_event_timestamps is None
            else max_event_timestamps
        ),
        min_amp=current.min_amp if min_amp is None else min_amp,
        min_imi=current.min_imi if min_imi is None else min_imi,
        nhood_radius=current.nhood_radius if nhood_radius is None else nhood_radius,
        x_up_threshold=(
            current.x_up_threshold if x_up_threshold is None else x_up_threshold
        ),
        x_down_threshold=(
            current.x_down_threshold if x_down_threshold is None else x_down_threshold
        ),
        dx_up_threshold=(
            current.dx_up_threshold if dx_up_threshold is None else dx_up_threshold
        ),
        dx_down_threshold=(
            current.dx_down_threshold
            if dx_down_threshold is None
            else dx_down_threshold
        ),
        eps_dx=current.eps_dx if eps_dx is None else eps_dx,
    )


def _resolve_variable_index(
    variable_names: Sequence[str],
    variable_name: str | None,
    *,
    default_index: int,
    parameter_name: str,
) -> int:
    if variable_name in (None, ""):
        return default_index

    names = tuple(variable_names)
    try:
        return names.index(variable_name)
    except ValueError as error:
        available_names = ", ".join(names) if names else "<none>"
        raise ValueError(
            f"Unknown {parameter_name} '{variable_name}'. Expected one of: {available_names}"
        ) from error


__all__ = ["Observer", "ObserverParams"]