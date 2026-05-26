from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..observers._definitions import ResolvedObserverSpec, resolve_observer_spec
from ..problem._core import ProblemInfo
from ..observers.types import EventOutputSettings, ObserverRuntimeSettings
from .models import Precision
from .runtime import OpenCLRuntime
from .structs import MatchedStruct, match_struct_dtype


@dataclass(frozen=True, slots=True)
class ObserverMetadata:
    observer_name: str
    feature_names: tuple[str, ...]
    observer_state_struct: MatchedStruct
    uses_two_pass: bool

    @property
    def n_features(self) -> int:
        return len(self.feature_names)

    @property
    def observer_state_nbytes(self) -> int:
        return self.observer_state_struct.dtype.itemsize


def get_observer_metadata(
    runtime: OpenCLRuntime,
    problem_info: ProblemInfo,
    observer_name: str,
    observer_runtime_settings: ObserverRuntimeSettings,
    event_output_settings: EventOutputSettings,
    precision: Precision,
    n_store_events: int | None = None,
    resolved_observer_spec: ResolvedObserverSpec | None = None,
) -> ObserverMetadata:
    resolved_spec = resolved_observer_spec
    if resolved_spec is None:
        resolved_n_store_events = (
            event_output_settings.max_event_timestamps
            if n_store_events is None
            else int(n_store_events)
        )
        resolved_spec = resolve_observer_spec(
            problem_info,
            observer_name,
            observer_runtime_settings,
            real_dtype=_real_dtype_for_precision(precision),
            event_output_settings=event_output_settings,
            n_store_events=resolved_n_store_events,
        )
    return ObserverMetadata(
        observer_name=resolved_spec.observer_name,
        feature_names=resolved_spec.feature_names,
        observer_state_struct=match_struct_dtype(
            runtime,
            resolved_spec.observer_state_struct_name,
            resolved_spec.observer_state_dtype,
        ),
        uses_two_pass=resolved_spec.uses_two_pass,
    )


def _real_dtype_for_precision(precision: Precision) -> np.dtype:
    return np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)