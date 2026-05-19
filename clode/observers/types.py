from __future__ import annotations

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
class _ObserverRuntimeSettings:
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


@dataclass(frozen=True, slots=True)
class _EventOutputSettings:
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
    """Configuration for built-in observer feature detection.

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
    def runtime_settings(self) -> _ObserverRuntimeSettings:
        return _ObserverRuntimeSettings(
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
    def event_output_settings(self) -> _EventOutputSettings:
        return _EventOutputSettings(
            max_event_timestamps=self.max_event_timestamps,
        )


__all__ = ["Observer", "ObserverParams"]