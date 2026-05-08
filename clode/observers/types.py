from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class Observer(Enum):
    basic = "basic"
    basic_all_variables = "basicall"
    local_max = "localmax"
    neighbourhood_1 = "nhood1"
    neighbourhood_2 = "nhood2"
    threshold_2 = "thresh2"


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


__all__ = ["Observer", "ObserverParams"]