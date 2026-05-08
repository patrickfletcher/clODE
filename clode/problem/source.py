from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path


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
    "RhsSource",
    "compute_rhs_digest",
    "create_rhs_source",
    "load_rhs_source",
]