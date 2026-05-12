from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


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


__all__ = ["Observer", "ObserverParams"]