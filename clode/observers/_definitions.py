from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

import numpy as np

from ..problem._core import ProblemInfo
from .types import Observer, ObserverParams


LayoutField = tuple[object, ...]
FeatureNameFactory = Callable[[ProblemInfo, ObserverParams], tuple[str, ...]]
LayoutFactory = Callable[[ProblemInfo, np.dtype, int], "ResolvedObserverLayout"]


@dataclass(frozen=True, slots=True)
class ResolvedObserverLayout:
    persistent_fields: tuple[LayoutField, ...]
    event_output_fields: tuple[LayoutField, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "persistent_fields",
            tuple(field for field in self.persistent_fields if field),
        )
        object.__setattr__(
            self,
            "event_output_fields",
            tuple(field for field in self.event_output_fields if field),
        )

    @property
    def all_fields(self) -> tuple[LayoutField, ...]:
        return self.persistent_fields + self.event_output_fields

    @property
    def dtype(self) -> np.dtype:
        return np.dtype(list(self.all_fields), align=True)


@dataclass(frozen=True, slots=True)
class ResolvedObserverSpec:
    definition: "ObserverDefinition"
    feature_names: tuple[str, ...]
    layout: ResolvedObserverLayout
    precision_name: str
    n_var: int
    n_aux: int
    n_store_events: int

    @property
    def n_features(self) -> int:
        return len(self.feature_names)

    @property
    def observer_name(self) -> str:
        return self.definition.observer_name

    @property
    def observer_data_dtype(self) -> np.dtype:
        return self.layout.dtype

    @property
    def observer_data_nbytes(self) -> int:
        return self.observer_data_dtype.itemsize

    @property
    def observer_data_struct_name(self) -> str:
        return (
            f"clode_observer_data_{self.observer_name}_{self.precision_name}"
            f"_v{self.n_var}_a{self.n_aux}_e{self.n_store_events}"
        )

    @property
    def uses_two_pass(self) -> bool:
        return self.definition.uses_two_pass

    @property
    def build_define(self) -> str:
        return self.definition.build_define

    @property
    def build_signature(self) -> tuple[str, int]:
        return (self.build_define, self.n_store_events)

    @property
    def layout_signature(self) -> tuple[str, int]:
        return (self.observer_data_struct_name, self.observer_data_nbytes)


@dataclass(frozen=True, slots=True)
class ObserverDefinition:
    observer_name: str
    build_define: str
    uses_two_pass: bool
    feature_name_factory: FeatureNameFactory
    layout_factory: LayoutFactory

    def get_feature_names(
        self,
        problem_info: ProblemInfo,
        observer_params: ObserverParams,
    ) -> tuple[str, ...]:
        return self.feature_name_factory(problem_info, observer_params)

    def resolve(
        self,
        problem_info: ProblemInfo,
        observer_params: ObserverParams,
        *,
        real_dtype: np.dtype,
        n_store_events: int | None = None,
    ) -> ResolvedObserverSpec:
        resolved_real_dtype = np.dtype(real_dtype)
        resolved_n_store_events = (
            observer_params.event_output_settings.max_event_timestamps
            if n_store_events is None
            else int(n_store_events)
        )
        return ResolvedObserverSpec(
            definition=self,
            feature_names=self.get_feature_names(problem_info, observer_params),
            layout=self.layout_factory(
                problem_info,
                resolved_real_dtype,
                resolved_n_store_events,
            ),
            precision_name=_precision_name(resolved_real_dtype),
            n_var=problem_info.num_var,
            n_aux=problem_info.num_aux,
            n_store_events=resolved_n_store_events,
        )


def get_observer_definition(observer: Observer | str) -> ObserverDefinition:
    observer_name = _normalize_observer_name(observer)
    try:
        return _OBSERVER_DEFINITIONS[observer_name]
    except KeyError as error:
        raise ValueError(f"Unsupported observer metadata request: {observer_name}") from error


def resolve_observer_spec(
    problem_info: ProblemInfo,
    observer: Observer | str,
    observer_params: ObserverParams,
    *,
    real_dtype: np.dtype,
    n_store_events: int | None = None,
) -> ResolvedObserverSpec:
    definition = get_observer_definition(observer)
    return definition.resolve(
        problem_info,
        observer_params,
        real_dtype=real_dtype,
        n_store_events=n_store_events,
    )


def _basic_feature_names(
    problem_info: ProblemInfo,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
    feature_var = _name_at(list(problem_info.vars), observer_params.f_var_ix)
    return (
        f"max {feature_var}",
        f"min {feature_var}",
        f"mean {feature_var}",
        f"max d{feature_var}/dt",
        f"min d{feature_var}/dt",
        "step count",
    )


def _basicall_feature_names(
    problem_info: ProblemInfo,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
    del observer_params
    names: list[str] = []
    for var_name in problem_info.vars:
        names.extend(
            [
                f"max {var_name}",
                f"min {var_name}",
                f"mean {var_name}",
                f"max d{var_name}/dt",
                f"min d{var_name}/dt",
            ]
        )
    for aux_name in problem_info.aux:
        names.extend(
            [
                f"max {aux_name}",
                f"min {aux_name}",
                f"mean {aux_name}",
            ]
        )
    names.append("step count")
    return tuple(names)


def _localmax_feature_names(
    problem_info: ProblemInfo,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
    names = [
        "max IMI",
        "min IMI",
        "mean IMI",
        "max amplitude",
        "min amplitude",
        "mean amplitude",
    ]
    for var_name in problem_info.vars:
        names.extend(
            [
                f"max {var_name}",
                f"min {var_name}",
                f"mean {var_name}",
                f"max d{var_name}/dt",
                f"min d{var_name}/dt",
            ]
        )
    for aux_name in problem_info.aux:
        names.extend(
            [
                f"max {aux_name}",
                f"min {aux_name}",
                f"mean {aux_name}",
            ]
        )
    for event_idx in range(observer_params.max_event_timestamps):
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


def _nhood1_feature_names(
    problem_info: ProblemInfo,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
    del observer_params
    names = [
        "max period",
        "min period",
        "mean period",
        "max peaks",
        "min peaks",
        "mean peaks",
    ]
    for var_name in problem_info.vars:
        names.extend(
            [
                f"max {var_name}",
                f"min {var_name}",
                f"mean {var_name}",
                f"max d{var_name}/dt",
                f"min d{var_name}/dt",
            ]
        )
    for aux_name in problem_info.aux:
        names.extend(
            [
                f"max {aux_name}",
                f"min {aux_name}",
                f"mean {aux_name}",
            ]
        )
    names.extend(["period count", "step count", "max dt", "min dt", "mean dt"])
    return tuple(names)


def _nhood2_feature_names(
    problem_info: ProblemInfo,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
    names = [
        "max period",
        "min period",
        "mean period",
        "max peaks",
        "min peaks",
        "mean peaks",
    ]
    for var_name in problem_info.vars:
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
    for aux_name in problem_info.aux:
        names.extend(
            [
                f"max {aux_name}",
                f"min {aux_name}",
                f"mean {aux_name}",
            ]
        )
    for event_idx in range(observer_params.max_event_timestamps):
        names.append(f"nhood event time {event_idx}")
    names.extend(["event count", "step count", "max dt", "min dt", "mean dt"])
    return tuple(names)


def _thresh2_feature_names(
    problem_info: ProblemInfo,
    observer_params: ObserverParams,
) -> tuple[str, ...]:
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
    for var_name in problem_info.vars:
        names.extend(
            [
                f"max {var_name}",
                f"min {var_name}",
                f"mean {var_name}",
                f"max d{var_name}/dt",
                f"min d{var_name}/dt",
            ]
        )
    for aux_name in problem_info.aux:
        names.extend(
            [
                f"max {aux_name}",
                f"min {aux_name}",
                f"mean {aux_name}",
            ]
        )
    for event_idx in range(observer_params.max_event_timestamps):
        names.extend(
            [
                f"up event time {event_idx}",
                f"down event time {event_idx}",
            ]
        )
    names.extend(["event count", "step count", "max dt", "min dt", "mean dt"])
    return tuple(names)


def _basic_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del n_store_events
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("xTrajectoryMax", real_dtype),
            _real_field("xTrajectoryMin", real_dtype),
            _real_field("xTrajectoryMean", real_dtype),
            _real_field("dxTrajectoryMax", real_dtype),
            _real_field("dxTrajectoryMin", real_dtype),
            _real_field("t_last", real_dtype),
            _real_field("t_start", real_dtype),
            _uint_field("stepcount"),
        )
    )


def _basicall_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del n_store_events
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("xTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMean", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("auxTrajectoryMax", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMin", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMean", real_dtype, problem_info.num_aux),
            _real_field("t_last", real_dtype),
            _real_field("t_start", real_dtype),
            _uint_field("stepcount"),
        )
    )


def _localmax_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    event_output_fields: list[LayoutField] = []
    if n_store_events > 0:
        event_output_fields.extend(
            [
                _real_field("tMaxList", real_dtype, n_store_events),
                _real_field("xMaxList", real_dtype, n_store_events),
                _real_field("tMinList", real_dtype, n_store_events),
                _real_field("xMinList", real_dtype, n_store_events),
            ]
        )
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("xbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("dxbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("xTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMean", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("auxTrajectoryMax", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMin", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMean", real_dtype, problem_info.num_aux),
            _real_field("IMI", real_dtype, 3),
            _real_field("amp", real_dtype, 3),
            _real_field("t_start", real_dtype),
            _real_field("tLastMax", real_dtype),
            _real_field("tLastMin", real_dtype),
            _real_field("xLastMin", real_dtype),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
        ),
        event_output_fields=tuple(event_output_fields),
    )


def _nhood1_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del n_store_events
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("xbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("dxbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("x0", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMean", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("auxTrajectoryMax", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMin", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMean", real_dtype, problem_info.num_aux),
            _real_field("nMaxima", real_dtype, 3),
            _real_field("period", real_dtype, 3),
            _real_field("stepDt", real_dtype, 3),
            _real_field("t_start", real_dtype),
            _real_field("tLastEvent", real_dtype),
            _real_field("thisNormXdiff", real_dtype),
            _real_field("lastNormXdiff", real_dtype),
            _uint_field("thisNMaxima"),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
            _uint_field("foundX0"),
            _uint_field("isInNhood"),
        )
    )


def _nhood2_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    event_output_fields: tuple[LayoutField, ...] = ()
    if n_store_events > 0:
        event_output_fields = (
            _real_field("tExitNhood", real_dtype, n_store_events),
        )
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("xbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("dxbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("x0", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMean", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryRange", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("auxTrajectoryMax", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMin", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMean", real_dtype, problem_info.num_aux),
            _real_field("nMaxima", real_dtype, 3),
            _real_field("period", real_dtype, 3),
            _real_field("stepDt", real_dtype, 3),
            _real_field("t_start", real_dtype),
            _real_field("tLastEvent", real_dtype),
            _real_field("xThreshold", real_dtype),
            _uint_field("thisNMaxima"),
            _uint_field("foundX0"),
            _uint_field("isInNhood"),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
        ),
        event_output_fields=event_output_fields,
    )


def _thresh2_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    event_output_fields: list[LayoutField] = []
    if n_store_events > 0:
        event_output_fields.extend(
            [
                _real_field("tUpTransition", real_dtype, n_store_events),
                _real_field("tDownTransition", real_dtype, n_store_events),
            ]
        )
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("xbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("dxbuffer", real_dtype, 3 * problem_info.num_var),
            _real_field("xTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMean", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("dxTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("auxTrajectoryMax", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMin", real_dtype, problem_info.num_aux),
            _real_field("auxTrajectoryMean", real_dtype, problem_info.num_aux),
            _real_field("nMaxima", real_dtype, 3),
            _real_field("period", real_dtype, 3),
            _real_field("upDuration", real_dtype, 3),
            _real_field("downDuration", real_dtype, 3),
            _real_field("duty", real_dtype, 3),
            _real_field("activeDip", real_dtype, 3),
            _real_field("stepDt", real_dtype, 3),
            _real_field("fVarUpstateMean", real_dtype),
            _real_field("fVarDownstateMean", real_dtype),
            _real_field("xGlobalMax", real_dtype),
            _real_field("xGlobalMin", real_dtype),
            _real_field("dxGlobalMax", real_dtype),
            _real_field("dxGlobalMin", real_dtype),
            _real_field("xUp", real_dtype),
            _real_field("xDown", real_dtype),
            _real_field("dxUp", real_dtype),
            _real_field("dxDown", real_dtype),
            _real_field("t_start", real_dtype),
            _real_field("tLastEvent", real_dtype),
            _real_field("tThisDown", real_dtype),
            _real_field("tLastMax", real_dtype),
            _real_field("tLastMin", real_dtype),
            _real_field("xLastMin", real_dtype),
            _uint_field("thisNMaxima"),
            _uint_field("stepcount"),
            _uint_field("eventcount"),
            _uint_field("inUpstate"),
        ),
        event_output_fields=tuple(event_output_fields),
    )


def _normalize_observer_name(observer: Observer | str) -> str:
    if isinstance(observer, Observer):
        return observer.value
    return str(observer)


def _name_at(names: list[str], index: int) -> str:
    if not names:
        raise ValueError("Observer feature naming requires at least one variable")
    clamped_index = min(max(int(index), 0), len(names) - 1)
    return names[clamped_index]


def _precision_name(real_dtype: np.dtype) -> str:
    if np.dtype(real_dtype) == np.dtype(np.float32):
        return "float"
    if np.dtype(real_dtype) == np.dtype(np.float64):
        return "double"
    raise ValueError(f"Unsupported observer precision dtype: {real_dtype}")


def _real_field(name: str, real_dtype: np.dtype, count: int = 1) -> LayoutField:
    if count <= 0:
        return ()
    if count == 1:
        return (name, real_dtype)
    return (name, real_dtype, (count,))


def _uint_field(name: str, count: int = 1) -> LayoutField:
    if count <= 0:
        return ()
    if count == 1:
        return (name, np.uint32)
    return (name, np.uint32, (count,))


_OBSERVER_DEFINITIONS = {
    "basic": ObserverDefinition(
        observer_name="basic",
        build_define="USE_OBSERVER_BASIC",
        uses_two_pass=False,
        feature_name_factory=_basic_feature_names,
        layout_factory=_basic_layout,
    ),
    "basicall": ObserverDefinition(
        observer_name="basicall",
        build_define="USE_OBSERVER_BASIC_ALLVAR",
        uses_two_pass=False,
        feature_name_factory=_basicall_feature_names,
        layout_factory=_basicall_layout,
    ),
    "localmax": ObserverDefinition(
        observer_name="localmax",
        build_define="USE_OBSERVER_LOCAL_MAX",
        uses_two_pass=False,
        feature_name_factory=_localmax_feature_names,
        layout_factory=_localmax_layout,
    ),
    "nhood1": ObserverDefinition(
        observer_name="nhood1",
        build_define="USE_OBSERVER_NEIGHBORHOOD_1",
        uses_two_pass=False,
        feature_name_factory=_nhood1_feature_names,
        layout_factory=_nhood1_layout,
    ),
    "nhood2": ObserverDefinition(
        observer_name="nhood2",
        build_define="USE_OBSERVER_NEIGHBORHOOD_2",
        uses_two_pass=True,
        feature_name_factory=_nhood2_feature_names,
        layout_factory=_nhood2_layout,
    ),
    "thresh2": ObserverDefinition(
        observer_name="thresh2",
        build_define="USE_OBSERVER_THRESHOLD_2",
        uses_two_pass=True,
        feature_name_factory=_thresh2_feature_names,
        layout_factory=_thresh2_layout,
    ),
}


__all__ = [
    "ObserverDefinition",
    "ResolvedObserverLayout",
    "ResolvedObserverSpec",
    "get_observer_definition",
    "resolve_observer_spec",
]