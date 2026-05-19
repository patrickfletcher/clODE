from __future__ import annotations

from dataclasses import dataclass
from enum import Enum


class StepperMethodKind(str, Enum):
    EXPLICIT = "explicit"
    IMPLICIT = "implicit"


class StepperStepSizeKind(str, Enum):
    FIXED = "fixed"
    ADAPTIVE = "adaptive"


class StepperNoiseKind(str, Enum):
    DETERMINISTIC = "deterministic"
    STOCHASTIC = "stochastic"


@dataclass(frozen=True, slots=True)
class StepperDefinition:
    stepper_name: str
    build_define: str
    method_kind: StepperMethodKind
    step_size_kind: StepperStepSizeKind
    noise_kind: StepperNoiseKind

    def __post_init__(self) -> None:
        if not self.stepper_name:
            raise ValueError("stepper_name must not be empty")
        if not self.build_define:
            raise ValueError("build_define must not be empty")

    @property
    def is_adaptive(self) -> bool:
        return self.step_size_kind is StepperStepSizeKind.ADAPTIVE

    @property
    def is_explicit(self) -> bool:
        return self.method_kind is StepperMethodKind.EXPLICIT

    @property
    def is_implicit(self) -> bool:
        return self.method_kind is StepperMethodKind.IMPLICIT

    @property
    def is_stochastic(self) -> bool:
        return self.noise_kind is StepperNoiseKind.STOCHASTIC


_STEPPER_DEFINITIONS = {
    "euler": StepperDefinition(
        stepper_name="euler",
        build_define="EXPLICIT_EULER",
        method_kind=StepperMethodKind.EXPLICIT,
        step_size_kind=StepperStepSizeKind.FIXED,
        noise_kind=StepperNoiseKind.DETERMINISTIC,
    ),
    "heun": StepperDefinition(
        stepper_name="heun",
        build_define="EXPLICIT_HEUN",
        method_kind=StepperMethodKind.EXPLICIT,
        step_size_kind=StepperStepSizeKind.FIXED,
        noise_kind=StepperNoiseKind.DETERMINISTIC,
    ),
    "rk4": StepperDefinition(
        stepper_name="rk4",
        build_define="EXPLICIT_RK4",
        method_kind=StepperMethodKind.EXPLICIT,
        step_size_kind=StepperStepSizeKind.FIXED,
        noise_kind=StepperNoiseKind.DETERMINISTIC,
    ),
    "bs23": StepperDefinition(
        stepper_name="bs23",
        build_define="EXPLICIT_BS23",
        method_kind=StepperMethodKind.EXPLICIT,
        step_size_kind=StepperStepSizeKind.ADAPTIVE,
        noise_kind=StepperNoiseKind.DETERMINISTIC,
    ),
    "dopri5": StepperDefinition(
        stepper_name="dopri5",
        build_define="EXPLICIT_DOPRI5",
        method_kind=StepperMethodKind.EXPLICIT,
        step_size_kind=StepperStepSizeKind.ADAPTIVE,
        noise_kind=StepperNoiseKind.DETERMINISTIC,
    ),
    "seuler": StepperDefinition(
        stepper_name="seuler",
        build_define="STOCHASTIC_EULER",
        method_kind=StepperMethodKind.EXPLICIT,
        step_size_kind=StepperStepSizeKind.FIXED,
        noise_kind=StepperNoiseKind.STOCHASTIC,
    ),
}


def get_stepper_definition(stepper_name: str) -> StepperDefinition:
    try:
        return _STEPPER_DEFINITIONS[stepper_name]
    except KeyError as error:
        raise ValueError(f"Unknown stepper: {stepper_name!r}") from error


def get_stepper_names() -> tuple[str, ...]:
    return tuple(_STEPPER_DEFINITIONS.keys())


def get_stepper_definitions() -> tuple[StepperDefinition, ...]:
    return tuple(_STEPPER_DEFINITIONS.values())
