from __future__ import annotations

from pathlib import Path

from ..observers._definitions import get_observer_definition
from .errors import UnsupportedObserverError, UnsupportedStepperError
from .models import KernelKind


class KernelRegistry:
    _STEPPERS = {
        "euler": "EXPLICIT_EULER",
        "heun": "EXPLICIT_HEUN",
        "rk4": "EXPLICIT_RK4",
        "bs23": "EXPLICIT_BS23",
        "dopri5": "EXPLICIT_DOPRI5",
        "seuler": "STOCHASTIC_EULER",
    }

    _ENTRYPOINTS = {
        KernelKind.TRANSIENT: ("transient.cl",),
        KernelKind.TRAJECTORY: ("transient.cl", "trajectory.cl"),
        KernelKind.FEATURES: (
            "transient.cl",
            "initializeObserver.cl",
            "features.cl",
        ),
    }

    _KERNEL_NAMES = {
        KernelKind.TRANSIENT: ("transient",),
        KernelKind.TRAJECTORY: ("transient", "trajectory"),
        KernelKind.FEATURES: ("transient", "initializeObserver", "features"),
    }

    def __init__(self, kernel_root: Path) -> None:
        self._kernel_root = kernel_root

    def get_stepper_define(self, stepper_name: str) -> str:
        self.validate_stepper(stepper_name)
        return self._STEPPERS[stepper_name]

    def get_observer_define(self, observer_name: str) -> str:
        return self._resolve_observer_definition(observer_name).build_define

    def get_entrypoint_paths(self, kernel_kind: KernelKind) -> tuple[Path, ...]:
        return tuple(self._kernel_root / entrypoint for entrypoint in self._ENTRYPOINTS[kernel_kind])

    def get_kernel_names(self, kernel_kind: KernelKind) -> tuple[str, ...]:
        return self._KERNEL_NAMES[kernel_kind]

    def validate_stepper(self, stepper_name: str) -> None:
        if stepper_name not in self._STEPPERS:
            raise UnsupportedStepperError(stepper_name)

    def validate_observer(self, observer_name: str | None) -> None:
        if observer_name is None:
            return
        self._resolve_observer_definition(observer_name)

    def _resolve_observer_definition(self, observer_name: str):
        try:
            return get_observer_definition(observer_name)
        except ValueError as error:
            raise UnsupportedObserverError(observer_name) from error
