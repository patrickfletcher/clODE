from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
import hashlib
from typing import Callable

import numpy as np

from ..problem._core import ProblemInfo
from .types import (
    EventOutputSettings,
    Observer,
    ObserverRuntimeSettings,
    SummaryObserverSelection,
)


LayoutField = tuple[object, ...]
FeatureNameFactory = Callable[
    [ProblemInfo, ObserverRuntimeSettings, int],
    tuple[str, ...],
]
LayoutFactory = Callable[[ProblemInfo, np.dtype, int], "ResolvedObserverLayout"]

_SUMMARY_BUILD_DEFINE = "USE_OBSERVER_SUMMARY"
_SUMMARY_KIND_ORDER: tuple[tuple[str, str], ...] = (
    ("state", "max"),
    ("state", "min"),
    ("state", "mean"),
    ("slope", "max"),
    ("slope", "min"),
    ("slope", "mean"),
    ("aux", "max"),
    ("aux", "min"),
    ("aux", "mean"),
)
_SUMMARY_LAYOUT_TAGS = {
    ("state", "max"): "sx",
    ("state", "min"): "sn",
    ("state", "mean"): "se",
    ("slope", "max"): "dx",
    ("slope", "min"): "dn",
    ("slope", "mean"): "de",
    ("aux", "max"): "ax",
    ("aux", "min"): "an",
    ("aux", "mean"): "ae",
}
_SUMMARY_COUNT_MACROS = {
    ("state", "max"): "CLODE_SUMMARY_STATE_MAX_COUNT",
    ("state", "min"): "CLODE_SUMMARY_STATE_MIN_COUNT",
    ("state", "mean"): "CLODE_SUMMARY_STATE_MEAN_COUNT",
    ("slope", "max"): "CLODE_SUMMARY_SLOPE_MAX_COUNT",
    ("slope", "min"): "CLODE_SUMMARY_SLOPE_MIN_COUNT",
    ("slope", "mean"): "CLODE_SUMMARY_SLOPE_MEAN_COUNT",
    ("aux", "max"): "CLODE_SUMMARY_AUX_MAX_COUNT",
    ("aux", "min"): "CLODE_SUMMARY_AUX_MIN_COUNT",
    ("aux", "mean"): "CLODE_SUMMARY_AUX_MEAN_COUNT",
}
_SUMMARY_INDEX_ARRAY_NAMES = {
    ("state", "max"): "CLODE_SUMMARY_STATE_MAX_INDICES",
    ("state", "min"): "CLODE_SUMMARY_STATE_MIN_INDICES",
    ("state", "mean"): "CLODE_SUMMARY_STATE_MEAN_INDICES",
    ("slope", "max"): "CLODE_SUMMARY_SLOPE_MAX_INDICES",
    ("slope", "min"): "CLODE_SUMMARY_SLOPE_MIN_INDICES",
    ("slope", "mean"): "CLODE_SUMMARY_SLOPE_MEAN_INDICES",
    ("aux", "max"): "CLODE_SUMMARY_AUX_MAX_INDICES",
    ("aux", "min"): "CLODE_SUMMARY_AUX_MIN_INDICES",
    ("aux", "mean"): "CLODE_SUMMARY_AUX_MEAN_INDICES",
}
_SUMMARY_OUTPUT_KIND_IDS = {
    kind: output_kind
    for output_kind, kind in enumerate(_SUMMARY_KIND_ORDER, start=1)
}


@dataclass(frozen=True, slots=True)
class ResolvedSummaryEntry:
    name: str
    index: int
    reductions: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ResolvedSummarySelection:
    state: tuple[ResolvedSummaryEntry, ...]
    aux: tuple[ResolvedSummaryEntry, ...]
    slope: tuple[ResolvedSummaryEntry, ...]
    build_variant: str
    source_preamble: str
    layout_name_tag: str

    def entries_for(self, group: str, reduction: str) -> tuple[ResolvedSummaryEntry, ...]:
        return tuple(
            entry for entry in getattr(self, group) if reduction in entry.reductions
        )


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
    build_variant: str = ""
    source_preamble: str = ""
    layout_name_tag: str | None = None

    @property
    def n_features(self) -> int:
        return len(self.feature_names)

    @property
    def observer_name(self) -> str:
        return self.definition.observer_name

    @property
    def observer_state_dtype(self) -> np.dtype:
        return self.layout.dtype

    @property
    def observer_state_nbytes(self) -> int:
        return self.observer_state_dtype.itemsize

    @property
    def observer_state_struct_name(self) -> str:
        name_tag = self.layout_name_tag or self.observer_name
        return (
            f"clode_observer_state_{name_tag}_{self.precision_name}"
            f"_v{self.n_var}_a{self.n_aux}_e{self.n_store_events}"
        )

    @property
    def uses_two_pass(self) -> bool:
        return self.definition.uses_two_pass

    @property
    def build_define(self) -> str:
        return self.definition.build_define

    @property
    def build_signature(self) -> tuple[str, int, str]:
        return (self.build_define, self.n_store_events, self.build_variant)

    @property
    def layout_signature(self) -> tuple[str, int]:
        return (self.observer_state_struct_name, self.observer_state_nbytes)


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
        observer_runtime_settings: ObserverRuntimeSettings,
        n_store_events: int,
        summary_selection: SummaryObserverSelection | None = None,
    ) -> tuple[str, ...]:
        if self.build_define == _SUMMARY_BUILD_DEFINE:
            return _summary_feature_names(
                _resolve_summary_selection(
                    self.observer_name,
                    problem_info,
                    observer_runtime_settings,
                    summary_selection,
                )
            )
        return self.feature_name_factory(
            problem_info,
            observer_runtime_settings,
            int(n_store_events),
        )

    def resolve(
        self,
        problem_info: ProblemInfo,
        observer_runtime_settings: ObserverRuntimeSettings,
        *,
        real_dtype: np.dtype,
        n_store_events: int,
        summary_selection: SummaryObserverSelection | None = None,
    ) -> ResolvedObserverSpec:
        resolved_real_dtype = np.dtype(real_dtype)
        resolved_n_store_events = int(n_store_events)
        resolved_summary_selection = (
            _resolve_summary_selection(
                self.observer_name,
                problem_info,
                observer_runtime_settings,
                summary_selection,
            )
            if self.build_define == _SUMMARY_BUILD_DEFINE
            else None
        )
        return ResolvedObserverSpec(
            definition=self,
            feature_names=(
                _summary_feature_names(resolved_summary_selection)
                if resolved_summary_selection is not None
                else self.get_feature_names(
                    problem_info,
                    observer_runtime_settings,
                    resolved_n_store_events,
                )
            ),
            layout=(
                _summary_layout(resolved_real_dtype, resolved_summary_selection)
                if resolved_summary_selection is not None
                else self.layout_factory(
                    problem_info,
                    resolved_real_dtype,
                    resolved_n_store_events,
                )
            ),
            precision_name=_precision_name(resolved_real_dtype),
            n_var=problem_info.num_var,
            n_aux=problem_info.num_aux,
            n_store_events=resolved_n_store_events,
            build_variant=(
                resolved_summary_selection.build_variant
                if resolved_summary_selection is not None
                else self.observer_name
            ),
            source_preamble=(
                resolved_summary_selection.source_preamble
                if resolved_summary_selection is not None
                else ""
            ),
            layout_name_tag=(
                resolved_summary_selection.layout_name_tag
                if resolved_summary_selection is not None
                else None
            ),
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
    observer_runtime_settings: ObserverRuntimeSettings,
    *,
    real_dtype: np.dtype,
    event_output_settings: EventOutputSettings,
    n_store_events: int | None = None,
    summary_selection: SummaryObserverSelection | None = None,
) -> ResolvedObserverSpec:
    definition = get_observer_definition(observer)
    if summary_selection is not None and definition.build_define != _SUMMARY_BUILD_DEFINE:
        raise ValueError(
            "summary_selection is only supported for summary observers"
        )
    resolved_n_store_events = (
        event_output_settings.max_event_timestamps
        if n_store_events is None
        else int(n_store_events)
    )
    return definition.resolve(
        problem_info,
        observer_runtime_settings,
        real_dtype=real_dtype,
        n_store_events=resolved_n_store_events,
        summary_selection=summary_selection,
    )


def _summary_feature_names(
    resolved_selection: ResolvedSummarySelection,
) -> tuple[str, ...]:
    feature_names: list[str] = []
    state_lookup = {entry.name: entry for entry in resolved_selection.state}
    slope_lookup = {entry.name: entry for entry in resolved_selection.slope}
    ordered_state_names = _ordered_unique(
        entry.name for entry in resolved_selection.state + resolved_selection.slope
    )

    for variable_name in ordered_state_names:
        state_entry = state_lookup.get(variable_name)
        slope_entry = slope_lookup.get(variable_name)
        if state_entry is not None:
            for reduction_name in ("max", "min", "mean"):
                if reduction_name in state_entry.reductions:
                    feature_names.append(f"{reduction_name} {variable_name}")
        if slope_entry is not None:
            for reduction_name in ("max", "min", "mean"):
                if reduction_name in slope_entry.reductions:
                    feature_names.append(f"{reduction_name} d{variable_name}/dt")

    for aux_entry in resolved_selection.aux:
        for reduction_name in ("max", "min", "mean"):
            if reduction_name in aux_entry.reductions:
                feature_names.append(f"{reduction_name} {aux_entry.name}")
    return tuple(feature_names)


def _localmax_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del observer_runtime_settings
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
    for event_idx in range(n_store_events):
        names.extend(
            [
                f"localmax event time {event_idx}",
                f"localmax event evar {event_idx}",
                f"localmin event time {event_idx}",
                f"localmin event evar {event_idx}",
            ]
        )
    names.append("event count")
    return tuple(names)


def _local_extremum_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del problem_info
    del observer_runtime_settings
    names: list[str] = []
    for event_idx in range(n_store_events):
        names.extend(
            [
                f"local_extremum event time {event_idx}",
                f"local_extremum event value {event_idx}",
            ]
        )
    names.append("event count")
    return tuple(names)


def _nhood1_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del observer_runtime_settings
    del n_store_events
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
    names.append("period count")
    return tuple(names)


def _nhood2_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del observer_runtime_settings
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
    for event_idx in range(n_store_events):
        names.append(f"nhood event time {event_idx}")
    names.append("event count")
    return tuple(names)


def _neighborhood_return_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del problem_info
    del observer_runtime_settings
    names = [
        f"neighborhood_return event time {event_idx}"
        for event_idx in range(n_store_events)
    ]
    names.append("event count")
    return tuple(names)


def _threshold_crossing_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del problem_info
    del observer_runtime_settings
    names = [f"threshold event time {event_idx}" for event_idx in range(n_store_events)]
    names.append("event count")
    return tuple(names)


def _normalized_threshold_crossing_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del problem_info
    del observer_runtime_settings
    names = [f"threshold event time {event_idx}" for event_idx in range(n_store_events)]
    names.append("event count")
    return tuple(names)


def _schmitt_trigger_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del problem_info
    del observer_runtime_settings
    names: list[str] = []
    for event_idx in range(n_store_events):
        names.extend(
            [
                f"up event time {event_idx}",
                f"down event time {event_idx}",
            ]
        )
    names.append("event count")
    return tuple(names)


def _threshold_2_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del observer_runtime_settings
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
    for event_idx in range(n_store_events):
        names.extend(
            [
                f"up event time {event_idx}",
                f"down event time {event_idx}",
            ]
        )
    names.append("event count")
    return tuple(names)


def _summary_layout(
    real_dtype: np.dtype,
    resolved_selection: ResolvedSummarySelection,
) -> ResolvedObserverLayout:
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field(
                "stateMax",
                real_dtype,
                len(resolved_selection.entries_for("state", "max")),
            ),
            _real_field(
                "stateMin",
                real_dtype,
                len(resolved_selection.entries_for("state", "min")),
            ),
            _real_field(
                "stateIntegral",
                real_dtype,
                len(resolved_selection.entries_for("state", "mean")),
            ),
            _real_field(
                "stateIntegralCorrection",
                real_dtype,
                len(resolved_selection.entries_for("state", "mean")),
            ),
            _real_field(
                "slopeMax",
                real_dtype,
                len(resolved_selection.entries_for("slope", "max")),
            ),
            _real_field(
                "slopeMin",
                real_dtype,
                len(resolved_selection.entries_for("slope", "min")),
            ),
            _real_field(
                "slopeIntegral",
                real_dtype,
                len(resolved_selection.entries_for("slope", "mean")),
            ),
            _real_field(
                "slopeIntegralCorrection",
                real_dtype,
                len(resolved_selection.entries_for("slope", "mean")),
            ),
            _real_field(
                "auxMax",
                real_dtype,
                len(resolved_selection.entries_for("aux", "max")),
            ),
            _real_field(
                "auxMin",
                real_dtype,
                len(resolved_selection.entries_for("aux", "min")),
            ),
            _real_field(
                "auxIntegral",
                real_dtype,
                len(resolved_selection.entries_for("aux", "mean")),
            ),
            _real_field(
                "auxIntegralCorrection",
                real_dtype,
                len(resolved_selection.entries_for("aux", "mean")),
            ),
            _real_field("elapsed_total", real_dtype),
            _real_field("elapsed_total_correction", real_dtype),
        )
    )


def _summary_placeholder_feature_names(
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    n_store_events: int,
) -> tuple[str, ...]:
    del problem_info
    del observer_runtime_settings
    del n_store_events
    raise RuntimeError("summary observers resolve feature names through selection metadata")


def _summary_placeholder_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del problem_info
    del real_dtype
    del n_store_events
    raise RuntimeError("summary observers resolve layouts through selection metadata")


def _resolve_summary_selection(
    observer_name: str,
    problem_info: ProblemInfo,
    observer_runtime_settings: ObserverRuntimeSettings,
    summary_selection: SummaryObserverSelection | None,
) -> ResolvedSummarySelection:
    if summary_selection is None:
        if observer_name == "basic":
            summary_selection = SummaryObserverSelection.single_variable(
                _name_at(list(problem_info.vars), observer_runtime_settings.f_var_ix)
            )
        else:
            summary_selection = SummaryObserverSelection.all_variables(
                problem_info.vars,
                problem_info.aux,
            )

    if summary_selection.is_empty:
        raise ValueError("Summary observer selection must include at least one reduction")

    resolved = ResolvedSummarySelection(
        state=_resolve_summary_group(problem_info.vars, summary_selection.state, "state"),
        aux=_resolve_summary_group(problem_info.aux, summary_selection.aux, "aux"),
        slope=_resolve_summary_group(problem_info.vars, summary_selection.slope, "slope"),
        build_variant="",
        source_preamble="",
        layout_name_tag="",
    )
    source_preamble = _summary_source_preamble(resolved)
    return ResolvedSummarySelection(
        state=resolved.state,
        aux=resolved.aux,
        slope=resolved.slope,
        build_variant=hashlib.sha256(source_preamble.encode("utf-8")).hexdigest(),
        source_preamble=source_preamble,
        layout_name_tag=_summary_layout_name_tag(resolved),
    )


def _resolve_summary_group(
    available_names: Iterable[str],
    selection_group: tuple[tuple[str, tuple[str, ...]], ...],
    group_name: str,
) -> tuple[ResolvedSummaryEntry, ...]:
    names = tuple(available_names)
    resolved_entries: list[ResolvedSummaryEntry] = []
    for variable_name, reductions in selection_group:
        try:
            variable_index = names.index(variable_name)
        except ValueError as error:
            available = ", ".join(names) if names else "<none>"
            raise ValueError(
                f"Unknown {group_name} summary variable '{variable_name}'. "
                f"Expected one of: {available}"
            ) from error
        resolved_entries.append(
            ResolvedSummaryEntry(variable_name, variable_index, tuple(reductions))
        )
    return tuple(resolved_entries)


def _summary_source_preamble(resolved_selection: ResolvedSummarySelection) -> str:
    lines: list[str] = []
    slot_lookups: dict[tuple[str, str], dict[int, int]] = {}

    for kind in _SUMMARY_KIND_ORDER:
        entries = resolved_selection.entries_for(*kind)
        lines.append(f"#define {_SUMMARY_COUNT_MACROS[kind]} {len(entries)}")
        if entries:
            slot_lookups[kind] = {
                entry.index: slot for slot, entry in enumerate(entries)
            }
            index_values = ", ".join(str(entry.index) for entry in entries)
            lines.append(
                f"__constant uint {_SUMMARY_INDEX_ARRAY_NAMES[kind]}[{len(entries)}] = {{{index_values}}};"
            )

    output_entries = _summary_output_entries(resolved_selection, slot_lookups)
    lines.append(f"#define CLODE_SUMMARY_OUTPUT_COUNT {len(output_entries)}")
    kind_values = ", ".join(
        str(_SUMMARY_OUTPUT_KIND_IDS[kind]) for kind, _slot in output_entries
    )
    slot_values = ", ".join(str(slot) for _kind, slot in output_entries)
    lines.append(
        f"__constant uint CLODE_SUMMARY_OUTPUT_KINDS[{len(output_entries)}] = {{{kind_values}}};"
    )
    lines.append(
        f"__constant uint CLODE_SUMMARY_OUTPUT_SLOTS[{len(output_entries)}] = {{{slot_values}}};"
    )
    return "\n".join(lines) + "\n"


def _summary_output_entries(
    resolved_selection: ResolvedSummarySelection,
    slot_lookups: dict[tuple[str, str], dict[int, int]],
) -> tuple[tuple[tuple[str, str], int], ...]:
    output_entries: list[tuple[tuple[str, str], int]] = []
    state_lookup = {entry.name: entry for entry in resolved_selection.state}
    slope_lookup = {entry.name: entry for entry in resolved_selection.slope}
    ordered_state_names = _ordered_unique(
        entry.name for entry in resolved_selection.state + resolved_selection.slope
    )

    for variable_name in ordered_state_names:
        state_entry = state_lookup.get(variable_name)
        slope_entry = slope_lookup.get(variable_name)
        if state_entry is not None:
            for reduction_name in ("max", "min", "mean"):
                if reduction_name in state_entry.reductions:
                    kind = ("state", reduction_name)
                    output_entries.append(
                        (kind, slot_lookups[kind][state_entry.index])
                    )
        if slope_entry is not None:
            for reduction_name in ("max", "min", "mean"):
                if reduction_name in slope_entry.reductions:
                    kind = ("slope", reduction_name)
                    output_entries.append(
                        (kind, slot_lookups[kind][slope_entry.index])
                    )

    for aux_entry in resolved_selection.aux:
        for reduction_name in ("max", "min", "mean"):
            if reduction_name in aux_entry.reductions:
                kind = ("aux", reduction_name)
                output_entries.append((kind, slot_lookups[kind][aux_entry.index]))

    return tuple(output_entries)


def _summary_layout_name_tag(
    resolved_selection: ResolvedSummarySelection,
) -> str:
    parts = ["summary"]
    for kind in _SUMMARY_KIND_ORDER:
        parts.append(
            f"{_SUMMARY_LAYOUT_TAGS[kind]}{len(resolved_selection.entries_for(*kind))}"
        )
    return "_".join(parts)


def _ordered_unique(values: Iterable[str]) -> tuple[str, ...]:
    ordered_values: list[str] = []
    seen: set[str] = set()
    for value in values:
        if value not in seen:
            ordered_values.append(value)
            seen.add(value)
    return tuple(ordered_values)


def _localmax_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("elapsedbuffer", real_dtype, 3),
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
            _real_field("elapsedTotal", real_dtype),
            _real_field("tLastMax", real_dtype),
            _real_field("elapsedLastMax", real_dtype),
            _real_field("tLastMin", real_dtype),
            _real_field("xLastMin", real_dtype),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
        ),
        event_output_fields=(
            _real_field("tMaxList", real_dtype, n_store_events),
            _real_field("xMaxList", real_dtype, n_store_events),
            _real_field("tMinList", real_dtype, n_store_events),
            _real_field("xMinList", real_dtype, n_store_events),
        ),
    )


def _local_extremum_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del problem_info
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("xbuffer", real_dtype, 3),
            _real_field("dxbuffer", real_dtype, 3),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
        ),
        event_output_fields=(
            _real_field("tEventList", real_dtype, n_store_events),
            _real_field("xEventList", real_dtype, n_store_events),
        ),
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
            _real_field("elapsedbuffer", real_dtype, 3),
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
            _real_field("elapsedTotal", real_dtype),
            _real_field("tLastEvent", real_dtype),
            _real_field("elapsedLastEvent", real_dtype),
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
            _real_field("tExitNhood", real_dtype, n_store_events),
            _real_field("nMaxima", real_dtype, 3),
            _real_field("period", real_dtype, 3),
            _real_field("elapsedTotal", real_dtype),
            _real_field("tLastEvent", real_dtype),
            _real_field("elapsedLastEvent", real_dtype),
            _real_field("xThreshold", real_dtype),
            _uint_field("thisNMaxima"),
            _uint_field("foundX0"),
            _uint_field("isInNhood"),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
        ),
    )


def _neighborhood_return_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("xbuffer", real_dtype, 2),
            _real_field("x0", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMax", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryMin", real_dtype, problem_info.num_var),
            _real_field("xTrajectoryRange", real_dtype, problem_info.num_var),
            _real_field("xThreshold", real_dtype),
            _uint_field("foundX0"),
            _uint_field("isInNhood"),
            _uint_field("eventcount"),
            _uint_field("stepcount"),
        ),
        event_output_fields=(
            _real_field("tEventList", real_dtype, n_store_events),
        ),
    )


def _threshold_crossing_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del problem_info
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 2),
            _real_field("xbuffer", real_dtype, 2),
            _real_field("xGlobalMax", real_dtype),
            _real_field("xGlobalMin", real_dtype),
            _uint_field("stepcount"),
            _uint_field("eventcount"),
        ),
        event_output_fields=(
            _real_field("tEventList", real_dtype, n_store_events),
        ),
    )


def _normalized_threshold_crossing_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del problem_info
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 2),
            _real_field("xbuffer", real_dtype, 2),
            _real_field("xGlobalMax", real_dtype),
            _real_field("xGlobalMin", real_dtype),
            _real_field("xThreshold", real_dtype),
            _uint_field("stepcount"),
            _uint_field("eventcount"),
        ),
        event_output_fields=(
            _real_field("tEventList", real_dtype, n_store_events),
        ),
    )


def _schmitt_trigger_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    del problem_info
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 2),
            _real_field("xbuffer", real_dtype, 2),
            _real_field("dxbuffer", real_dtype, 2),
            _real_field("xGlobalMax", real_dtype),
            _real_field("xGlobalMin", real_dtype),
            _real_field("dxGlobalMax", real_dtype),
            _real_field("dxGlobalMin", real_dtype),
            _real_field("xUp", real_dtype),
            _real_field("xDown", real_dtype),
            _real_field("dxUp", real_dtype),
            _real_field("dxDown", real_dtype),
            _uint_field("stepcount"),
            _uint_field("eventcount"),
            _uint_field("inUpstate"),
        ),
        event_output_fields=(
            _real_field("tUpTransition", real_dtype, n_store_events),
            _real_field("tDownTransition", real_dtype, n_store_events),
        ),
    )


def _threshold_2_layout(
    problem_info: ProblemInfo,
    real_dtype: np.dtype,
    n_store_events: int,
) -> ResolvedObserverLayout:
    return ResolvedObserverLayout(
        persistent_fields=(
            _real_field("tbuffer", real_dtype, 3),
            _real_field("elapsedbuffer", real_dtype, 3),
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
            _real_field("tUpTransition", real_dtype, n_store_events),
            _real_field("tDownTransition", real_dtype, n_store_events),
            _real_field("nMaxima", real_dtype, 3),
            _real_field("period", real_dtype, 3),
            _real_field("upDuration", real_dtype, 3),
            _real_field("downDuration", real_dtype, 3),
            _real_field("duty", real_dtype, 3),
            _real_field("activeDip", real_dtype, 3),
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
            _real_field("elapsedTotal", real_dtype),
            _real_field("tLastEvent", real_dtype),
            _real_field("elapsedLastEvent", real_dtype),
            _real_field("tThisDown", real_dtype),
            _real_field("elapsedThisDown", real_dtype),
            _real_field("tLastMax", real_dtype),
            _real_field("tLastMin", real_dtype),
            _real_field("xLastMin", real_dtype),
            _uint_field("thisNMaxima"),
            _uint_field("stepcount"),
            _uint_field("eventcount"),
            _uint_field("inUpstate"),
        ),
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
    "summary": ObserverDefinition(
        observer_name="summary",
        build_define=_SUMMARY_BUILD_DEFINE,
        uses_two_pass=False,
        feature_name_factory=_summary_placeholder_feature_names,
        layout_factory=_summary_placeholder_layout,
    ),
    "basic": ObserverDefinition(
        observer_name="basic",
        build_define=_SUMMARY_BUILD_DEFINE,
        uses_two_pass=False,
        feature_name_factory=_summary_placeholder_feature_names,
        layout_factory=_summary_placeholder_layout,
    ),
    "basicall": ObserverDefinition(
        observer_name="basicall",
        build_define=_SUMMARY_BUILD_DEFINE,
        uses_two_pass=False,
        feature_name_factory=_summary_placeholder_feature_names,
        layout_factory=_summary_placeholder_layout,
    ),
    "local_extremum": ObserverDefinition(
        observer_name="local_extremum",
        build_define="USE_OBSERVER_LOCAL_EXTREMUM",
        uses_two_pass=False,
        feature_name_factory=_local_extremum_feature_names,
        layout_factory=_local_extremum_layout,
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
    "neighborhood_return": ObserverDefinition(
        observer_name="neighborhood_return",
        build_define="USE_OBSERVER_NEIGHBORHOOD_RETURN",
        uses_two_pass=True,
        feature_name_factory=_neighborhood_return_feature_names,
        layout_factory=_neighborhood_return_layout,
    ),
    "threshold_crossing": ObserverDefinition(
        observer_name="threshold_crossing",
        build_define="USE_OBSERVER_THRESHOLD_CROSSING",
        uses_two_pass=False,
        feature_name_factory=_threshold_crossing_feature_names,
        layout_factory=_threshold_crossing_layout,
    ),
    "normalized_threshold_crossing": ObserverDefinition(
        observer_name="normalized_threshold_crossing",
        build_define="USE_OBSERVER_NORMALIZED_THRESHOLD_CROSSING",
        uses_two_pass=True,
        feature_name_factory=_normalized_threshold_crossing_feature_names,
        layout_factory=_normalized_threshold_crossing_layout,
    ),
    "schmitt_trigger": ObserverDefinition(
        observer_name="schmitt_trigger",
        build_define="USE_OBSERVER_SCHMITT_TRIGGER",
        uses_two_pass=False,
        feature_name_factory=_schmitt_trigger_feature_names,
        layout_factory=_schmitt_trigger_layout,
    ),
    "normalized_schmitt_trigger": ObserverDefinition(
        observer_name="normalized_schmitt_trigger",
        build_define="USE_OBSERVER_NORMALIZED_SCHMITT_TRIGGER",
        uses_two_pass=True,
        feature_name_factory=_schmitt_trigger_feature_names,
        layout_factory=_schmitt_trigger_layout,
    ),
    "threshold_2": ObserverDefinition(
        observer_name="threshold_2",
        build_define="USE_OBSERVER_THRESHOLD_2",
        uses_two_pass=True,
        feature_name_factory=_threshold_2_feature_names,
        layout_factory=_threshold_2_layout,
    ),
}


__all__ = [
    "ObserverDefinition",
    "ResolvedObserverLayout",
    "ResolvedObserverSpec",
    "get_observer_definition",
    "resolve_observer_spec",
]