from __future__ import annotations

from ..problem._core import ProblemInfo
from ._definitions import get_observer_definition
from .types import Observer, ObserverParams, SummaryObserverSelection


def is_two_pass_observer(observer: Observer | str) -> bool:
    return get_observer_definition(observer).uses_two_pass


def get_observer_feature_names(
    problem_info: ProblemInfo,
    observer: Observer | str,
    observer_params: ObserverParams,
    summary_selection: SummaryObserverSelection | None = None,
) -> tuple[str, ...]:
    return get_observer_definition(observer).get_feature_names(
        problem_info,
        observer_params.runtime_settings,
        observer_params.event_output_settings.max_event_timestamps,
        summary_selection=summary_selection,
    )


__all__ = ["get_observer_feature_names", "is_two_pass_observer"]