from __future__ import annotations

from ..problem.definition import ProblemInfo
from .types import Observer, ObserverParams


_TWO_PASS_OBSERVERS = frozenset({"nhood2", "thresh2"})


def is_two_pass_observer(observer: Observer | str) -> bool:
    return _normalize_observer_name(observer) in _TWO_PASS_OBSERVERS


def get_observer_feature_names(
    problem_info: ProblemInfo,
    observer: Observer | str,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
    observer_name = _normalize_observer_name(observer)
    var_names = list(problem_info.vars)
    aux_names = list(problem_info.aux)
    feature_var = _name_at(var_names, observer_params.f_var_ix)
    n_store_events = observer_params.max_event_timestamps

    if observer_name == "basic":
        return (
            f"max {feature_var}",
            f"min {feature_var}",
            f"mean {feature_var}",
            f"max d{feature_var}/dt",
            f"min d{feature_var}/dt",
            "step count",
        )

    if observer_name == "basicall":
        names: list[str] = []
        for var_name in var_names:
            names.extend(
                [
                    f"max {var_name}",
                    f"min {var_name}",
                    f"mean {var_name}",
                    f"max d{var_name}/dt",
                    f"min d{var_name}/dt",
                ]
            )
        for aux_name in aux_names:
            names.extend(
                [
                    f"max {aux_name}",
                    f"min {aux_name}",
                    f"mean {aux_name}",
                ]
            )
        names.append("step count")
        return tuple(names)

    if observer_name == "localmax":
        names = [
            "max IMI",
            "min IMI",
            "mean IMI",
            "max amplitude",
            "min amplitude",
            "mean amplitude",
        ]
        for var_name in var_names:
            names.extend(
                [
                    f"max {var_name}",
                    f"min {var_name}",
                    f"mean {var_name}",
                    f"max d{var_name}/dt",
                    f"min d{var_name}/dt",
                ]
            )
        for aux_name in aux_names:
            names.extend(
                [
                    f"max {aux_name}",
                    f"min {aux_name}",
                    f"mean {aux_name}",
                ]
            )
        for event_idx in range(n_store_events):
            names.extend(
                [
                    f"localmax event time {event_idx}",
                    f"localmax event evar {event_idx}",
                    f"localmin event time {event_idx}",
                    f"localmin event evar {event_idx}",
                ]
            )
        names.extend(["event count", "step count"])
        return tuple(names)

    if observer_name == "nhood1":
        names = [
            "max period",
            "min period",
            "mean period",
            "max peaks",
            "min peaks",
            "mean peaks",
        ]
        for var_name in var_names:
            names.extend(
                [
                    f"max {var_name}",
                    f"min {var_name}",
                    f"mean {var_name}",
                    f"max d{var_name}/dt",
                    f"min d{var_name}/dt",
                ]
            )
        for aux_name in aux_names:
            names.extend(
                [
                    f"max {aux_name}",
                    f"min {aux_name}",
                    f"mean {aux_name}",
                ]
            )
        names.extend(["period count", "step count", "max dt", "min dt", "mean dt"])
        return tuple(names)

    if observer_name == "nhood2":
        names = [
            "max period",
            "min period",
            "mean period",
            "max peaks",
            "min peaks",
            "mean peaks",
        ]
        for var_name in var_names:
            names.extend(
                [
                    f"max {var_name}",
                    f"min {var_name}",
                    f"mean {var_name}",
                    f"range {var_name}",
                    f"nhood center {var_name}",
                    f"max d{var_name}/dt",
                    f"min d{var_name}/dt",
                ]
            )
        for aux_name in aux_names:
            names.extend(
                [
                    f"max {aux_name}",
                    f"min {aux_name}",
                    f"mean {aux_name}",
                ]
            )
        for event_idx in range(n_store_events):
            names.append(f"nhood event time {event_idx}")
        names.extend(["event count", "step count", "max dt", "min dt", "mean dt"])
        return tuple(names)

    if observer_name == "thresh2":
        names = [
            "max period",
            "min period",
            "mean period",
            "max peaks",
            "min peaks",
            "mean peaks",
            "max upDuration",
            "min upDuration",
            "mean upDuration",
            "max downDuration",
            "min downDuration",
            "mean downDuration",
            "max duty",
            "min duty",
            "mean duty",
            "max activeDip",
            "min activeDip",
            "mean activeDip",
        ]
        for var_name in var_names:
            names.extend(
                [
                    f"max {var_name}",
                    f"min {var_name}",
                    f"mean {var_name}",
                    f"max d{var_name}/dt",
                    f"min d{var_name}/dt",
                ]
            )
        for aux_name in aux_names:
            names.extend(
                [
                    f"max {aux_name}",
                    f"min {aux_name}",
                    f"mean {aux_name}",
                ]
            )
        for event_idx in range(n_store_events):
            names.extend(
                [
                    f"up event time {event_idx}",
                    f"down event time {event_idx}",
                ]
            )
        names.extend(["event count", "step count", "max dt", "min dt", "mean dt"])
        return tuple(names)

    raise ValueError(f"Unsupported observer metadata request: {observer_name}")


def _normalize_observer_name(observer: Observer | str) -> str:
    if isinstance(observer, Observer):
        return observer.value
    return str(observer)


def _name_at(names: list[str], index: int) -> str:
    if not names:
        raise ValueError("Observer feature naming requires at least one variable")
    clamped_index = min(max(int(index), 0), len(names) - 1)
    return names[clamped_index]


__all__ = ["get_observer_feature_names", "is_two_pass_observer"]