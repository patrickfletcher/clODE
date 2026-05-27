from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

from ..problem._core import RhsSource

if TYPE_CHECKING:
    from ..problem._core import ProblemInfo


OPENCL_BACKEND_VERSION = "1"


class Precision(str, Enum):
    SINGLE = "single"
    DOUBLE = "double"

    @classmethod
    def from_single_precision(cls, single_precision: bool) -> Precision:
        return cls.SINGLE if single_precision else cls.DOUBLE


class KernelKind(str, Enum):
    TRANSIENT = "transient"
    TRAJECTORY = "trajectory"
    FEATURES = "features"


@dataclass(frozen=True, slots=True)
class ProblemShape:
    n_var: int
    n_par: int
    n_aux: int
    n_wiener: int

    def __post_init__(self) -> None:
        for field_name in ("n_var", "n_par", "n_aux", "n_wiener"):
            value = getattr(self, field_name)
            if value < 0:
                raise ValueError(f"{field_name} must be non-negative")

    @classmethod
    def from_problem_info(cls, problem_info: ProblemInfo) -> ProblemShape:
        return cls(
            n_var=problem_info.num_var,
            n_par=problem_info.num_par,
            n_aux=problem_info.num_aux,
            n_wiener=problem_info.num_noise,
        )


@dataclass(frozen=True, slots=True)
class BuildKey:
    backend_version: str
    kernel_kind: KernelKind
    precision: Precision
    stepper_name: str
    stepper_define: str
    observer_define: str | None
    observer_signature: str | None
    problem_shape: ProblemShape
    n_store_events: int
    rhs_digest: str
    kernel_tree_digest: str
    debug_build: bool = False

    def __post_init__(self) -> None:
        if not self.backend_version:
            raise ValueError("backend_version must not be empty")
        if not self.stepper_name:
            raise ValueError("stepper_name must not be empty")
        if not self.stepper_define:
            raise ValueError("stepper_define must not be empty")
        if self.observer_define == "":
            raise ValueError("observer_define must be None or a non-empty string")
        if self.observer_signature == "":
            raise ValueError("observer_signature must be None or a non-empty string")
        if self.n_store_events < 0:
            raise ValueError("n_store_events must be non-negative")
        if self.observer_define is None and self.n_store_events != 0:
            raise ValueError("n_store_events requires an observer_define")
        if self.observer_define is None and self.observer_signature is not None:
            raise ValueError("observer_signature requires an observer_define")
        if self.observer_define is not None and self.observer_signature is None:
            raise ValueError("observer_signature is required when observer_define is set")
        if not self.rhs_digest:
            raise ValueError("rhs_digest must not be empty")
        if not self.kernel_tree_digest:
            raise ValueError("kernel_tree_digest must not be empty")


@dataclass(frozen=True, slots=True)
class SourceBundle:
    build_key: BuildKey
    source_text: str
    build_options: tuple[str, ...]
    kernel_names: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.source_text:
            raise ValueError("source_text must not be empty")
        if not self.kernel_names:
            raise ValueError("kernel_names must not be empty")


@dataclass(slots=True)
class ProgramBundle:
    build_key: BuildKey
    source_bundle: SourceBundle
    program: object
    kernels: dict[str, object]

    def __post_init__(self) -> None:
        if self.build_key != self.source_bundle.build_key:
            raise ValueError("Program bundle build_key must match the source bundle build_key")
        missing_kernels = [
            kernel_name
            for kernel_name in self.source_bundle.kernel_names
            if kernel_name not in self.kernels
        ]
        if missing_kernels:
            raise ValueError(
                "Program bundle is missing kernel handles for: "
                + ", ".join(missing_kernels)
            )
