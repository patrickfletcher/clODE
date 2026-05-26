"""Compatibility barrel retained for legacy type-bundle imports."""

from __future__ import annotations

from .observers.types import ObserverParams
from .simulation.params import SolverParams

__all__ = [
    "ObserverParams",
    "SolverParams",
]