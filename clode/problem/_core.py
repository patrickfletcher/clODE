from __future__ import annotations

from dataclasses import dataclass, field
import hashlib
from pathlib import Path


@dataclass(slots=True)
class ProblemInfo:
    """Static shape metadata derived from one modeled ODE problem."""

    src_file: str = ""
    vars: list[str] = field(default_factory=list)
    pars: list[str] = field(default_factory=list)
    aux: list[str] = field(default_factory=list)
    num_noise: int = 1

    def __post_init__(self) -> None:
        self.vars = list(self.vars)
        self.pars = list(self.pars)
        self.aux = list(self.aux)
        self.num_noise = int(self.num_noise)
        if self.num_noise < 0:
            raise ValueError("num_noise must be non-negative")

    @property
    def num_var(self) -> int:
        return len(self.vars)

    @property
    def num_par(self) -> int:
        return len(self.pars)

    @property
    def num_aux(self) -> int:
        return len(self.aux)


@dataclass(frozen=True)
class RhsSource:
    origin_label: str
    text: str
    digest: str


def compute_rhs_digest(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def create_rhs_source(origin_label: str, text: str) -> RhsSource:
    return RhsSource(
        origin_label=origin_label,
        text=text,
        digest=compute_rhs_digest(text),
    )


def load_rhs_source(source_path: str) -> RhsSource:
    text = Path(source_path).read_text(encoding="utf-8")
    return create_rhs_source(source_path, text)


__all__ = [
    "ProblemInfo",
    "RhsSource",
    "compute_rhs_digest",
    "create_rhs_source",
    "load_rhs_source",
]