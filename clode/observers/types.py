from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from enum import Enum, IntEnum


SummaryReductionSpec = str | Sequence[str]
NormalizedSummaryGroup = tuple[tuple[str, tuple[str, ...]], ...]

_SUMMARY_REDUCTION_ORDER = ("max", "min", "mean")
_SUMMARY_REDUCTION_NAMES = frozenset(_SUMMARY_REDUCTION_ORDER)


class EventDirection(IntEnum):
    """Direction selector for threshold-style event crossings."""

    rising = 0
    falling = 1
    either = 2

# common observer settings
_DEFAULT_EVENT_VAR_INDEX = 0
_DEFAULT_FEATURE_VAR_INDEX = 0

# observers with event detection
_DEFAULT_MAX_EVENT_COUNT = 100
_DEFAULT_MAX_EVENT_TIMESTAMPS = 0

# oscillatory limiters for event acceptance
_DEFAULT_MIN_AMP = 0.0
_DEFAULT_MIN_IMI = 0.0  # not used? maybe in local extremum to limit events too close together? probably a better way to do it though

# threshold crossing observers
_DEFAULT_THRESHOLD = 0.0
_DEFAULT_EVENT_DIRECTION = EventDirection.rising

# neighborhood observers
_DEFAULT_NHOOD_RADIUS = 0.05

# Schmitt trigger observers
_DEFAULT_X_UP_THRESHOLD = 0.3
_DEFAULT_X_DOWN_THRESHOLD = 0.2
# Derivative thresholds remain part of ObserverParams for compatibility plumbing,
# but semantic Schmitt observers do not use them.
_DEFAULT_DX_UP_THRESHOLD = 0.0
_DEFAULT_DX_DOWN_THRESHOLD = 0.0

_DEFAULT_EPS_DX = 0.0  #unused?


@dataclass(frozen=True, slots=True)
class ObserverRuntimeSettings:
    """Internal value object for observer runtime thresholds and selector indices."""

    e_var_ix: int = _DEFAULT_EVENT_VAR_INDEX
    f_var_ix: int = _DEFAULT_FEATURE_VAR_INDEX
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    event_direction: EventDirection | str | int = _DEFAULT_EVENT_DIRECTION
    min_amp: float = _DEFAULT_MIN_AMP
    min_imi: float = _DEFAULT_MIN_IMI
    nhood_radius: float = _DEFAULT_NHOOD_RADIUS
    x_up_threshold: float = _DEFAULT_X_UP_THRESHOLD
    x_down_threshold: float = _DEFAULT_X_DOWN_THRESHOLD
    dx_up_threshold: float = _DEFAULT_DX_UP_THRESHOLD
    dx_down_threshold: float = _DEFAULT_DX_DOWN_THRESHOLD
    eps_dx: float = _DEFAULT_EPS_DX

    def __post_init__(self) -> None:
        object.__setattr__(self, "e_var_ix", int(self.e_var_ix))
        object.__setattr__(self, "f_var_ix", int(self.f_var_ix))
        object.__setattr__(self, "max_event_count", int(self.max_event_count))
        object.__setattr__(
            self,
            "event_direction",
            _normalize_event_direction(self.event_direction),
        )
        object.__setattr__(self, "min_amp", float(self.min_amp))
        object.__setattr__(self, "min_imi", float(self.min_imi))
        object.__setattr__(self, "nhood_radius", float(self.nhood_radius))
        object.__setattr__(self, "x_up_threshold", float(self.x_up_threshold))
        object.__setattr__(self, "x_down_threshold", float(self.x_down_threshold))
        object.__setattr__(self, "dx_up_threshold", float(self.dx_up_threshold))
        object.__setattr__(self, "dx_down_threshold", float(self.dx_down_threshold))
        object.__setattr__(self, "eps_dx", float(self.eps_dx))

    def to_observer_params(
        self, event_output_settings: EventOutputSettings
    ) -> ObserverParams:
        return ObserverParams(
            e_var_ix=self.e_var_ix,
            f_var_ix=self.f_var_ix,
            max_event_count=self.max_event_count,
            event_direction=self.event_direction,
            max_event_timestamps=event_output_settings.max_event_timestamps,
            min_amp=self.min_amp,
            min_imi=self.min_imi,
            nhood_radius=self.nhood_radius,
            x_up_threshold=self.x_up_threshold,
            x_down_threshold=self.x_down_threshold,
            dx_up_threshold=self.dx_up_threshold,
            dx_down_threshold=self.dx_down_threshold,
            eps_dx=self.eps_dx,
        )


@dataclass(frozen=True, slots=True)
class EventOutputSettings:
    """Internal value object for retained observer event-output policy."""

    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS

    def __post_init__(self) -> None:
        object.__setattr__(self, "max_event_timestamps", int(self.max_event_timestamps))


class Observer(Enum):
    """Built-in observer modes available to `FeatureSimulator`.

    Each entry names a distinct semantic observer family. There are no
    compatibility aliases — use the canonical names directly.
    """

    summary = "summary"
    local_max = "localmax"
    threshold_crossing = "threshold_crossing"
    normalized_threshold_crossing = "normalized_threshold_crossing"
    schmitt_trigger = "schmitt_trigger"
    normalized_schmitt_trigger = "normalized_schmitt_trigger"
    normalized_neighborhood_return = "normalized_neighborhood_return"


@dataclass(frozen=True, slots=True)
class ThresholdCrossingConfig:
    """Semantic config for absolute and fractional threshold-crossing observers."""

    event_var: str = ""
    feature_var: str = ""
    threshold: float = _DEFAULT_THRESHOLD
    direction: EventDirection | str | int = _DEFAULT_EVENT_DIRECTION
    min_amp: float = _DEFAULT_MIN_AMP
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "event_var",
            _normalize_optional_variable_name(self.event_var, parameter_name="event_var"),
        )
        object.__setattr__(
            self,
            "feature_var",
            _normalize_optional_variable_name(
                self.feature_var,
                parameter_name="feature_var",
            ),
        )
        object.__setattr__(self, "threshold", float(self.threshold))
        object.__setattr__(self, "direction", _normalize_event_direction(self.direction))
        object.__setattr__(self, "min_amp", float(self.min_amp))
        object.__setattr__(self, "max_event_count", int(self.max_event_count))
        object.__setattr__(self, "max_event_timestamps", int(self.max_event_timestamps))

    @property
    def supported_observers(self) -> tuple[Observer, Observer]:
        return (
            Observer.threshold_crossing,
            Observer.normalized_threshold_crossing,
        )

    def to_observer_params(self, variable_names: Sequence[str]) -> ObserverParams:
        threshold = float(self.threshold)
        return ObserverParams(
            e_var_ix=_resolve_variable_index(
                variable_names,
                self.event_var or None,
                default_index=_DEFAULT_EVENT_VAR_INDEX,
                parameter_name="event_var",
            ),
            f_var_ix=_resolve_variable_index(
                variable_names,
                self.feature_var or None,
                default_index=_DEFAULT_FEATURE_VAR_INDEX,
                parameter_name="feature_var",
            ),
            max_event_count=self.max_event_count,
            event_direction=self.direction,
            max_event_timestamps=self.max_event_timestamps,
            min_amp=self.min_amp,
            x_up_threshold=threshold,
            x_down_threshold=threshold,
        )

    @classmethod
    def from_observer_params(
        cls,
        variable_names: Sequence[str],
        observer_params: ObserverParams,
    ) -> ThresholdCrossingConfig:
        return cls(
            event_var=_resolve_variable_name(
                variable_names,
                observer_params.e_var_ix,
                parameter_name="e_var_ix",
            ),
            feature_var=_resolve_variable_name(
                variable_names,
                observer_params.f_var_ix,
                parameter_name="f_var_ix",
            ),
            threshold=observer_params.x_up_threshold,
            direction=observer_params.event_direction,
            min_amp=observer_params.min_amp,
            max_event_count=observer_params.max_event_count,
            max_event_timestamps=observer_params.max_event_timestamps,
        )


@dataclass(frozen=True, slots=True)
class SchmittTriggerConfig:
    """Semantic config for the lean absolute and normalized Schmitt families."""

    event_var: str = ""
    feature_var: str = ""
    x_up_threshold: float = _DEFAULT_X_UP_THRESHOLD
    x_down_threshold: float = _DEFAULT_X_DOWN_THRESHOLD
    min_amp: float = _DEFAULT_MIN_AMP
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "event_var",
            _normalize_optional_variable_name(self.event_var, parameter_name="event_var"),
        )
        object.__setattr__(
            self,
            "feature_var",
            _normalize_optional_variable_name(
                self.feature_var,
                parameter_name="feature_var",
            ),
        )
        object.__setattr__(self, "x_up_threshold", float(self.x_up_threshold))
        object.__setattr__(self, "x_down_threshold", float(self.x_down_threshold))
        object.__setattr__(self, "min_amp", float(self.min_amp))
        object.__setattr__(self, "max_event_count", int(self.max_event_count))
        object.__setattr__(self, "max_event_timestamps", int(self.max_event_timestamps))
        if self.x_up_threshold < self.x_down_threshold:
            raise ValueError(
                "SchmittTriggerConfig requires x_up_threshold >= x_down_threshold"
            )

    @property
    def supported_observers(self) -> tuple[Observer, Observer]:
        return (
            Observer.schmitt_trigger,
            Observer.normalized_schmitt_trigger,
        )

    def to_observer_params(self, variable_names: Sequence[str]) -> ObserverParams:
        return ObserverParams(
            e_var_ix=_resolve_variable_index(
                variable_names,
                self.event_var or None,
                default_index=_DEFAULT_EVENT_VAR_INDEX,
                parameter_name="event_var",
            ),
            f_var_ix=_resolve_variable_index(
                variable_names,
                self.feature_var or None,
                default_index=_DEFAULT_FEATURE_VAR_INDEX,
                parameter_name="feature_var",
            ),
            max_event_count=self.max_event_count,
            event_direction=EventDirection.either,
            max_event_timestamps=self.max_event_timestamps,
            min_amp=self.min_amp,
            x_up_threshold=self.x_up_threshold,
            x_down_threshold=self.x_down_threshold,
            dx_up_threshold=_DEFAULT_DX_UP_THRESHOLD,
            dx_down_threshold=_DEFAULT_DX_DOWN_THRESHOLD,
        )

    @classmethod
    def from_observer_params(
        cls,
        variable_names: Sequence[str],
        observer_params: ObserverParams,
    ) -> SchmittTriggerConfig:
        return cls(
            event_var=_resolve_variable_name(
                variable_names,
                observer_params.e_var_ix,
                parameter_name="e_var_ix",
            ),
            feature_var=_resolve_variable_name(
                variable_names,
                observer_params.f_var_ix,
                parameter_name="f_var_ix",
            ),
            x_up_threshold=observer_params.x_up_threshold,
            x_down_threshold=observer_params.x_down_threshold,
            min_amp=observer_params.min_amp,
            max_event_count=observer_params.max_event_count,
            max_event_timestamps=observer_params.max_event_timestamps,
        )


@dataclass(frozen=True, slots=True)
class LocalMaximumConfig:
    """Semantic config for the canonical local-maximum observer.

    This observer is maxima-only and uses one variable for both event triggering
    and feature measurement (eVarIx == fVarIx).
    """

    event_var: str = ""
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS
    min_amp: float = _DEFAULT_MIN_AMP

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "event_var",
            _normalize_optional_variable_name(self.event_var, parameter_name="event_var"),
        )
        object.__setattr__(self, "min_amp", float(self.min_amp))
        object.__setattr__(self, "max_event_count", int(self.max_event_count))
        object.__setattr__(self, "max_event_timestamps", int(self.max_event_timestamps))

    @property
    def supported_observers(self) -> tuple[Observer]:
        return (Observer.local_max,)

    def to_observer_params(self, variable_names: Sequence[str]) -> ObserverParams:
        variable_ix = _resolve_variable_index(
            variable_names,
            self.event_var or None,
            default_index=_DEFAULT_EVENT_VAR_INDEX,
            parameter_name="event_var",
        )
        return ObserverParams(
            e_var_ix=variable_ix,
            f_var_ix=variable_ix,
            min_amp=self.min_amp,
            max_event_count=self.max_event_count,
            max_event_timestamps=self.max_event_timestamps,
        )

    @classmethod
    def from_observer_params(
        cls,
        variable_names: Sequence[str],
        observer_params: ObserverParams,
    ) -> LocalMaximumConfig:
        return cls(
            event_var=_resolve_variable_name(
                variable_names,
                observer_params.e_var_ix,
                parameter_name="e_var_ix",
            ),
            min_amp=observer_params.min_amp,
            max_event_count=observer_params.max_event_count,
            max_event_timestamps=observer_params.max_event_timestamps,
        )


@dataclass(frozen=True, slots=True)
class NeighborhoodReturnConfig:
    """Semantic config for the lean normalized neighborhood-return observer."""

    event_var: str = ""
    feature_var: str = ""
    anchor_threshold: float = _DEFAULT_X_DOWN_THRESHOLD
    radius: float = _DEFAULT_NHOOD_RADIUS
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS
    min_amp: float = _DEFAULT_MIN_AMP

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "event_var",
            _normalize_optional_variable_name(self.event_var, parameter_name="event_var"),
        )
        object.__setattr__(
            self,
            "feature_var",
            _normalize_optional_variable_name(
                self.feature_var,
                parameter_name="feature_var",
            ),
        )
        object.__setattr__(self, "anchor_threshold", float(self.anchor_threshold))
        object.__setattr__(self, "radius", float(self.radius))
        object.__setattr__(self, "min_amp", float(self.min_amp))
        object.__setattr__(self, "max_event_count", int(self.max_event_count))
        object.__setattr__(self, "max_event_timestamps", int(self.max_event_timestamps))
        _require_unit_interval(
            self.anchor_threshold,
            parameter_name="anchor_threshold",
        )
        _require_positive_float(self.radius, parameter_name="radius")

    @property
    def supported_observers(self) -> tuple[Observer]:
        return (Observer.normalized_neighborhood_return,)

    def to_observer_params(self, variable_names: Sequence[str]) -> ObserverParams:
        return ObserverParams(
            e_var_ix=_resolve_variable_index(
                variable_names,
                self.event_var or None,
                default_index=_DEFAULT_EVENT_VAR_INDEX,
                parameter_name="event_var",
            ),
            f_var_ix=_resolve_variable_index(
                variable_names,
                self.feature_var or None,
                default_index=_DEFAULT_FEATURE_VAR_INDEX,
                parameter_name="feature_var",
            ),
            min_amp=self.min_amp,
            max_event_count=self.max_event_count,
            max_event_timestamps=self.max_event_timestamps,
            nhood_radius=self.radius,
            x_down_threshold=self.anchor_threshold,
        )

    @classmethod
    def from_observer_params(
        cls,
        variable_names: Sequence[str],
        observer_params: ObserverParams,
    ) -> NeighborhoodReturnConfig:
        return cls(
            event_var=_resolve_variable_name(
                variable_names,
                observer_params.e_var_ix,
                parameter_name="e_var_ix",
            ),
            feature_var=_resolve_variable_name(
                variable_names,
                observer_params.f_var_ix,
                parameter_name="f_var_ix",
            ),
            anchor_threshold=observer_params.x_down_threshold,
            radius=observer_params.nhood_radius,
            min_amp=observer_params.min_amp,
            max_event_count=observer_params.max_event_count,
            max_event_timestamps=observer_params.max_event_timestamps,
        )


ObserverConfiguration = (
    ThresholdCrossingConfig
    | SchmittTriggerConfig
    | LocalMaximumConfig
    | NeighborhoodReturnConfig
)


@dataclass(frozen=True, slots=True)
class SummaryObserverSelection:
    """Explicit variable and reduction selection for summary observers.

    The mapping for each group preserves insertion order. Values may be a single
    reduction name or a sequence drawn from `"max"`, `"min"`, and `"mean"`.
    """

    state: NormalizedSummaryGroup = ()
    aux: NormalizedSummaryGroup = ()
    slope: NormalizedSummaryGroup = ()

    def __init__(
        self,
        *,
        state: Mapping[str, SummaryReductionSpec] | None = None,
        aux: Mapping[str, SummaryReductionSpec] | None = None,
        slope: Mapping[str, SummaryReductionSpec] | None = None,
    ) -> None:
        object.__setattr__(self, "state", _normalize_summary_group(state, "state"))
        object.__setattr__(self, "aux", _normalize_summary_group(aux, "aux"))
        object.__setattr__(self, "slope", _normalize_summary_group(slope, "slope"))

    @property
    def is_empty(self) -> bool:
        return not (self.state or self.aux or self.slope)

    @classmethod
    def single_variable(cls, variable_name: str) -> SummaryObserverSelection:
        return cls(
            state={variable_name: ("max", "min", "mean")},
            slope={variable_name: ("max", "min")},
        )

    @classmethod
    def all_variables(
        cls,
        variable_names: Sequence[str],
        aux_names: Sequence[str] = (),
    ) -> SummaryObserverSelection:
        state = {
            variable_name: ("max", "min", "mean")
            for variable_name in variable_names
        }
        aux = {aux_name: ("max", "min", "mean") for aux_name in aux_names}
        slope = {variable_name: ("max", "min") for variable_name in variable_names}
        return cls(state=state, aux=aux, slope=slope)


@dataclass(slots=True)
class ObserverParams:
    """Public compatibility bundle for built-in observer configuration.

    This object is retained for the current user-facing API. Internal code should
    prefer `ObserverRuntimeSettings` and `EventOutputSettings` when the narrower
    owner model is sufficient.

    Attributes:
        e_var_ix: Index of the variable used for event detection.
        f_var_ix: Index of the variable used for feature readout.
        max_event_count: Maximum number of events to accumulate.
        event_direction: Crossing direction used by threshold-style observers.
        max_event_timestamps: Maximum number of event timestamps to retain.
        min_amp: Minimum amplitude threshold for event acceptance.
        min_imi: Minimum inter-event interval.
        nhood_radius: Neighborhood radius used by neighborhood observers.
        x_up_threshold: Rising threshold for threshold-style observers.
        x_down_threshold: Falling threshold for threshold-style observers.
        dx_up_threshold: Rising derivative threshold for threshold-style observers.
        dx_down_threshold: Falling derivative threshold for threshold-style observers.
        eps_dx: Small derivative tolerance used around threshold crossings.
    """

    e_var_ix: int = _DEFAULT_EVENT_VAR_INDEX
    f_var_ix: int = _DEFAULT_FEATURE_VAR_INDEX
    max_event_count: int = _DEFAULT_MAX_EVENT_COUNT
    event_direction: EventDirection | str | int = _DEFAULT_EVENT_DIRECTION
    max_event_timestamps: int = _DEFAULT_MAX_EVENT_TIMESTAMPS
    min_amp: float = _DEFAULT_MIN_AMP
    min_imi: float = _DEFAULT_MIN_IMI
    nhood_radius: float = _DEFAULT_NHOOD_RADIUS
    x_up_threshold: float = _DEFAULT_X_UP_THRESHOLD
    x_down_threshold: float = _DEFAULT_X_DOWN_THRESHOLD
    dx_up_threshold: float = _DEFAULT_DX_UP_THRESHOLD
    dx_down_threshold: float = _DEFAULT_DX_DOWN_THRESHOLD
    eps_dx: float = _DEFAULT_EPS_DX

    def __post_init__(self) -> None:
        self.e_var_ix = int(self.e_var_ix)
        self.f_var_ix = int(self.f_var_ix)
        self.max_event_count = int(self.max_event_count)
        self.event_direction = _normalize_event_direction(self.event_direction)
        self.max_event_timestamps = int(self.max_event_timestamps)
        self.min_amp = float(self.min_amp)
        self.min_imi = float(self.min_imi)
        self.nhood_radius = float(self.nhood_radius)
        self.x_up_threshold = float(self.x_up_threshold)
        self.x_down_threshold = float(self.x_down_threshold)
        self.dx_up_threshold = float(self.dx_up_threshold)
        self.dx_down_threshold = float(self.dx_down_threshold)
        self.eps_dx = float(self.eps_dx)

    @property
    def runtime_settings(self) -> ObserverRuntimeSettings:
        return ObserverRuntimeSettings(
            e_var_ix=self.e_var_ix,
            f_var_ix=self.f_var_ix,
            max_event_count=self.max_event_count,
            event_direction=self.event_direction,
            min_amp=self.min_amp,
            min_imi=self.min_imi,
            nhood_radius=self.nhood_radius,
            x_up_threshold=self.x_up_threshold,
            x_down_threshold=self.x_down_threshold,
            dx_up_threshold=self.dx_up_threshold,
            dx_down_threshold=self.dx_down_threshold,
            eps_dx=self.eps_dx,
        )

    @property
    def event_output_settings(self) -> EventOutputSettings:
        return EventOutputSettings(
            max_event_timestamps=self.max_event_timestamps,
        )


def _copy_observer_params(observer_params: ObserverParams) -> ObserverParams:
    return ObserverParams(
        e_var_ix=observer_params.e_var_ix,
        f_var_ix=observer_params.f_var_ix,
        max_event_count=observer_params.max_event_count,
        event_direction=observer_params.event_direction,
        max_event_timestamps=observer_params.max_event_timestamps,
        min_amp=observer_params.min_amp,
        min_imi=observer_params.min_imi,
        nhood_radius=observer_params.nhood_radius,
        x_up_threshold=observer_params.x_up_threshold,
        x_down_threshold=observer_params.x_down_threshold,
        dx_up_threshold=observer_params.dx_up_threshold,
        dx_down_threshold=observer_params.dx_down_threshold,
        eps_dx=observer_params.eps_dx,
    )


def _resolve_observer_params(
    variable_names: Sequence[str],
    *,
    observer_params: ObserverParams | None = None,
    base_params: ObserverParams | None = None,
    event_var: str | None = None,
    feature_var: str | None = None,
    max_event_count: int | None = None,
    event_direction: EventDirection | str | int | None = None,
    max_event_timestamps: int | None = None,
    min_amp: float | None = None,
    min_imi: float | None = None,
    nhood_radius: float | None = None,
    x_up_threshold: float | None = None,
    x_down_threshold: float | None = None,
    dx_up_threshold: float | None = None,
    dx_down_threshold: float | None = None,
    eps_dx: float | None = None,
) -> ObserverParams:
    if observer_params is not None:
        return _copy_observer_params(observer_params)

    current = ObserverParams() if base_params is None else _copy_observer_params(base_params)
    return ObserverParams(
        e_var_ix=_resolve_variable_index(
            variable_names,
            event_var,
            default_index=current.e_var_ix,
            parameter_name="event_var",
        ),
        f_var_ix=_resolve_variable_index(
            variable_names,
            feature_var,
            default_index=current.f_var_ix,
            parameter_name="feature_var",
        ),
        max_event_count=current.max_event_count if max_event_count is None else max_event_count,
        event_direction=(
            current.event_direction
            if event_direction is None
            else _normalize_event_direction(event_direction)
        ),
        max_event_timestamps=(
            current.max_event_timestamps
            if max_event_timestamps is None
            else max_event_timestamps
        ),
        min_amp=current.min_amp if min_amp is None else min_amp,
        min_imi=current.min_imi if min_imi is None else min_imi,
        nhood_radius=current.nhood_radius if nhood_radius is None else nhood_radius,
        x_up_threshold=(
            current.x_up_threshold if x_up_threshold is None else x_up_threshold
        ),
        x_down_threshold=(
            current.x_down_threshold if x_down_threshold is None else x_down_threshold
        ),
        dx_up_threshold=(
            current.dx_up_threshold if dx_up_threshold is None else dx_up_threshold
        ),
        dx_down_threshold=(
            current.dx_down_threshold
            if dx_down_threshold is None
            else dx_down_threshold
        ),
        eps_dx=current.eps_dx if eps_dx is None else eps_dx,
    )


def _resolve_observer_configuration(
    variable_names: Sequence[str],
    *,
    observer: Observer | None,
    observer_configuration: ObserverConfiguration,
) -> tuple[Observer, ObserverParams]:
    supported_observers = observer_configuration.supported_observers
    if observer is None:
        if len(supported_observers) == 1:
            resolved_observer = supported_observers[0]
            resolved_params = observer_configuration.to_observer_params(variable_names)
            _validate_observer_params_for_observer(resolved_observer, resolved_params)
            return resolved_observer, resolved_params
        raise ValueError(
            "observer_configuration requires an explicit observer because the "
            "same config surface is shared across multiple observer variants"
        )
    if observer not in supported_observers:
        supported_names = ", ".join(supported.name for supported in supported_observers)
        raise ValueError(
            "observer_configuration does not match the requested observer. "
            f"Expected one of: {supported_names}; got '{observer.name}'."
        )
    resolved_params = observer_configuration.to_observer_params(variable_names)
    _validate_observer_params_for_observer(observer, resolved_params)
    return observer, resolved_params


def _observer_configuration_from_params(
    variable_names: Sequence[str],
    *,
    observer: Observer,
    observer_params: ObserverParams,
) -> ObserverConfiguration | None:
    if observer in (
        Observer.threshold_crossing,
        Observer.normalized_threshold_crossing,
    ):
        return ThresholdCrossingConfig.from_observer_params(
            variable_names,
            observer_params,
        )
    if observer in (
        Observer.schmitt_trigger,
        Observer.normalized_schmitt_trigger,
    ):
        return SchmittTriggerConfig.from_observer_params(
            variable_names,
            observer_params,
        )
    if observer is Observer.local_max:
        return LocalMaximumConfig.from_observer_params(
            variable_names,
            observer_params,
        )
    if observer is Observer.normalized_neighborhood_return:
        return NeighborhoodReturnConfig.from_observer_params(
            variable_names,
            observer_params,
        )
    return None


def _resolve_variable_index(
    variable_names: Sequence[str],
    variable_name: str | None,
    *,
    default_index: int,
    parameter_name: str,
) -> int:
    if variable_name in (None, ""):
        return default_index

    names = tuple(variable_names)
    try:
        return names.index(variable_name)
    except ValueError as error:
        available_names = ", ".join(names) if names else "<none>"
        raise ValueError(
            f"Unknown {parameter_name} '{variable_name}'. Expected one of: {available_names}"
        ) from error


def _resolve_variable_name(
    variable_names: Sequence[str],
    index: int,
    *,
    parameter_name: str,
) -> str:
    names = tuple(variable_names)
    if not names:
        return ""
    if not 0 <= int(index) < len(names):
        raise ValueError(
            f"{parameter_name} index {index} is out of range for variables: {', '.join(names)}"
        )
    return names[int(index)]


def _normalize_optional_variable_name(
    variable_name: str | None,
    *,
    parameter_name: str,
) -> str:
    if variable_name in (None, ""):
        return ""
    if not isinstance(variable_name, str):
        raise ValueError(f"{parameter_name} must be a string when provided")
    normalized_name = variable_name.strip()
    if not normalized_name:
        return ""
    return normalized_name


def _normalize_event_direction(
    event_direction: EventDirection | str | int,
) -> EventDirection:
    if isinstance(event_direction, EventDirection):
        return event_direction

    if isinstance(event_direction, str):
        normalized_direction = event_direction.strip().lower()
        aliases = {
            "up": "rising",
            "increasing": "rising",
            "down": "falling",
            "decreasing": "falling",
            "both": "either",
            "any": "either",
        }
        normalized_direction = aliases.get(normalized_direction, normalized_direction)
        try:
            return EventDirection[normalized_direction]
        except KeyError as error:
            allowed = ", ".join(direction.name for direction in EventDirection)
            raise ValueError(
                f"Unsupported event_direction '{event_direction}'. Expected one of: {allowed}"
            ) from error

    try:
        return EventDirection(int(event_direction))
    except (TypeError, ValueError) as error:
        allowed = ", ".join(direction.name for direction in EventDirection)
        raise ValueError(
            f"Unsupported event_direction '{event_direction}'. Expected one of: {allowed}"
        ) from error


def _validate_observer_params_for_observer(
    observer: Observer,
    observer_params: ObserverParams,
) -> None:
    if observer is Observer.normalized_threshold_crossing:
        _require_unit_interval(
            observer_params.x_up_threshold,
            parameter_name="threshold",
        )
        return

    if observer in (Observer.schmitt_trigger, Observer.normalized_schmitt_trigger):
        if (
            observer_params.dx_up_threshold != _DEFAULT_DX_UP_THRESHOLD
            or observer_params.dx_down_threshold != _DEFAULT_DX_DOWN_THRESHOLD
        ):
            raise ValueError(
                "Semantic Schmitt observers do not support dx thresholds. "
                "Use value thresholds only for semantic Schmitt observers."
            )
        if observer_params.x_up_threshold < observer_params.x_down_threshold:
            raise ValueError(
                "Schmitt observers require x_up_threshold >= x_down_threshold"
            )
        if observer is Observer.normalized_schmitt_trigger:
            _require_unit_interval(
                observer_params.x_up_threshold,
                parameter_name="x_up_threshold",
            )
            _require_unit_interval(
                observer_params.x_down_threshold,
                parameter_name="x_down_threshold",
            )
        return

    if observer is Observer.normalized_neighborhood_return:
        _require_unit_interval(
            observer_params.x_down_threshold,
            parameter_name="anchor_threshold",
        )
        _require_positive_float(
            observer_params.nhood_radius,
            parameter_name="radius",
        )


def _require_unit_interval(value: float, *, parameter_name: str) -> None:
    if not 0.0 <= float(value) <= 1.0:
        raise ValueError(f"{parameter_name} must be between 0.0 and 1.0 inclusive")


def _require_positive_float(value: float, *, parameter_name: str) -> None:
    if float(value) <= 0.0:
        raise ValueError(f"{parameter_name} must be greater than 0.0")


def _normalize_summary_group(
    group: Mapping[str, SummaryReductionSpec] | None,
    group_name: str,
) -> NormalizedSummaryGroup:
    if group is None:
        return ()

    normalized: list[tuple[str, tuple[str, ...]]] = []
    for variable_name, reduction_spec in group.items():
        if not isinstance(variable_name, str) or not variable_name:
            raise ValueError(f"{group_name} selection keys must be non-empty strings")
        reductions = _normalize_reduction_spec(
            reduction_spec,
            group_name=group_name,
            variable_name=variable_name,
        )
        normalized.append((variable_name, reductions))
    return tuple(normalized)


def _normalize_reduction_spec(
    reduction_spec: SummaryReductionSpec,
    *,
    group_name: str,
    variable_name: str,
) -> tuple[str, ...]:
    raw_reductions = (
        (reduction_spec,)
        if isinstance(reduction_spec, str)
        else tuple(reduction_spec)
    )
    if not raw_reductions:
        raise ValueError(
            f"{group_name} selection for '{variable_name}' must include at least one reduction"
        )

    normalized: list[str] = []
    seen: set[str] = set()
    for reduction_name in raw_reductions:
        normalized_name = str(reduction_name).strip().lower()
        if normalized_name not in _SUMMARY_REDUCTION_NAMES:
            allowed = ", ".join(_SUMMARY_REDUCTION_ORDER)
            raise ValueError(
                f"Unsupported {group_name} reduction '{reduction_name}' for '{variable_name}'. "
                f"Expected one of: {allowed}"
            )
        if normalized_name not in seen:
            normalized.append(normalized_name)
            seen.add(normalized_name)
    return tuple(normalized)


__all__ = [
    "EventDirection",
    "LocalMaximumConfig",
    "NeighborhoodReturnConfig",
    "Observer",
    "ObserverParams",
    "SchmittTriggerConfig",
    "SummaryObserverSelection",
    "ThresholdCrossingConfig",
]