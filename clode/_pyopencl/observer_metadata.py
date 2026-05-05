from __future__ import annotations

from dataclasses import dataclass

from clode.cpp.clode_cpp_wrapper import ObserverParams, ProblemInfo
import numpy as np

from .models import Precision, ProblemShape
from .runtime import OpenCLRuntime
from .structs import MatchedStruct, match_struct_dtype


_TWO_PASS_OBSERVERS = frozenset({"nhood2", "thresh2"})


@dataclass(frozen=True, slots=True)
class ObserverMetadata:
    observer_name: str
    feature_names: tuple[str, ...]
    observer_data_struct: MatchedStruct
    uses_two_pass: bool

    @property
    def n_features(self) -> int:
        return len(self.feature_names)

    @property
    def observer_data_nbytes(self) -> int:
        return self.observer_data_struct.dtype.itemsize


def get_observer_metadata(
    runtime: OpenCLRuntime,
    problem_info: ProblemInfo,
    observer_name: str,
    observer_params: ObserverParams,
    precision: Precision,
) -> ObserverMetadata:
    shape = ProblemShape.from_problem_info(problem_info)
    n_store_events = observer_params.max_event_timestamps
    base_dtype = _observer_data_base_dtype(
        observer_name,
        shape,
        precision,
        n_store_events,
    )
    struct_name = _observer_data_struct_name(
        observer_name,
        shape,
        precision,
        n_store_events,
    )
    return ObserverMetadata(
        observer_name=observer_name,
        feature_names=tuple(
            _feature_names_for_observer(problem_info, observer_name, observer_params)
        ),
        observer_data_struct=match_struct_dtype(runtime, struct_name, base_dtype),
        uses_two_pass=observer_name in _TWO_PASS_OBSERVERS,
    )


def is_two_pass_observer(observer_name: str) -> bool:
    return observer_name in _TWO_PASS_OBSERVERS


def _feature_names_for_observer(
    problem_info: ProblemInfo,
    observer_name: str,
    observer_params: ObserverParams,
) -> list[str]:
    var_names = list(problem_info.vars)
    aux_names = list(problem_info.aux)
    feature_var = _name_at(var_names, observer_params.f_var_ix)
    n_store_events = observer_params.max_event_timestamps

    if observer_name == "basic":
        return [
            f"max {feature_var}",
            f"min {feature_var}",
            f"mean {feature_var}",
            f"max d{feature_var}/dt",
            f"min d{feature_var}/dt",
            "step count",
        ]

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
        return names

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
        return names

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
        return names

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
        return names

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
        return names

    raise ValueError(f"Unsupported observer metadata request: {observer_name}")


def _observer_data_base_dtype(
    observer_name: str,
    shape: ProblemShape,
    precision: Precision,
    n_store_events: int,
) -> np.dtype:
    real_dtype = np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)
    fields: list[tuple[object, ...]] = []

    if observer_name == "basic":
        _append_real_field(fields, "xTrajectoryMax", real_dtype)
        _append_real_field(fields, "xTrajectoryMin", real_dtype)
        _append_real_field(fields, "xTrajectoryMean", real_dtype)
        _append_real_field(fields, "dxTrajectoryMax", real_dtype)
        _append_real_field(fields, "dxTrajectoryMin", real_dtype)
        _append_real_field(fields, "t_last", real_dtype)
        _append_real_field(fields, "t_start", real_dtype)
        _append_uint_field(fields, "stepcount")
        return np.dtype(fields, align=True)

    if observer_name == "basicall":
        _append_real_field(fields, "xTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMean", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "auxTrajectoryMax", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMin", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMean", real_dtype, shape.n_aux)
        _append_real_field(fields, "t_last", real_dtype)
        _append_real_field(fields, "t_start", real_dtype)
        _append_uint_field(fields, "stepcount")
        return np.dtype(fields, align=True)

    if observer_name == "localmax":
        _append_real_field(fields, "tbuffer", real_dtype, 3)
        _append_real_field(fields, "xbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "dxbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "xTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMean", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "auxTrajectoryMax", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMin", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMean", real_dtype, shape.n_aux)
        _append_real_field(fields, "tMaxList", real_dtype, n_store_events)
        _append_real_field(fields, "xMaxList", real_dtype, n_store_events)
        _append_real_field(fields, "tMinList", real_dtype, n_store_events)
        _append_real_field(fields, "xMinList", real_dtype, n_store_events)
        _append_real_field(fields, "IMI", real_dtype, 3)
        _append_real_field(fields, "amp", real_dtype, 3)
        _append_real_field(fields, "t_start", real_dtype)
        _append_real_field(fields, "tLastMax", real_dtype)
        _append_real_field(fields, "tLastMin", real_dtype)
        _append_real_field(fields, "xLastMin", real_dtype)
        _append_uint_field(fields, "eventcount")
        _append_uint_field(fields, "stepcount")
        return np.dtype(fields, align=True)

    if observer_name == "nhood1":
        _append_real_field(fields, "tbuffer", real_dtype, 3)
        _append_real_field(fields, "xbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "dxbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "x0", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMean", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "auxTrajectoryMax", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMin", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMean", real_dtype, shape.n_aux)
        _append_real_field(fields, "nMaxima", real_dtype, 3)
        _append_real_field(fields, "period", real_dtype, 3)
        _append_real_field(fields, "stepDt", real_dtype, 3)
        _append_real_field(fields, "t_start", real_dtype)
        _append_real_field(fields, "tLastEvent", real_dtype)
        _append_real_field(fields, "thisNormXdiff", real_dtype)
        _append_real_field(fields, "lastNormXdiff", real_dtype)
        _append_uint_field(fields, "thisNMaxima")
        _append_uint_field(fields, "eventcount")
        _append_uint_field(fields, "stepcount")
        _append_uint_field(fields, "foundX0")
        _append_uint_field(fields, "isInNhood")
        return np.dtype(fields, align=True)

    if observer_name == "nhood2":
        _append_real_field(fields, "tbuffer", real_dtype, 3)
        _append_real_field(fields, "xbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "dxbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "x0", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMean", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryRange", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "auxTrajectoryMax", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMin", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMean", real_dtype, shape.n_aux)
        _append_real_field(fields, "tExitNhood", real_dtype, n_store_events)
        _append_real_field(fields, "nMaxima", real_dtype, 3)
        _append_real_field(fields, "period", real_dtype, 3)
        _append_real_field(fields, "stepDt", real_dtype, 3)
        _append_real_field(fields, "t_start", real_dtype)
        _append_real_field(fields, "tLastEvent", real_dtype)
        _append_real_field(fields, "xThreshold", real_dtype)
        _append_uint_field(fields, "thisNMaxima")
        _append_uint_field(fields, "foundX0")
        _append_uint_field(fields, "isInNhood")
        _append_uint_field(fields, "eventcount")
        _append_uint_field(fields, "stepcount")
        return np.dtype(fields, align=True)

    if observer_name == "thresh2":
        _append_real_field(fields, "tbuffer", real_dtype, 3)
        _append_real_field(fields, "xbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "dxbuffer", real_dtype, 3 * shape.n_var)
        _append_real_field(fields, "xTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "xTrajectoryMean", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMax", real_dtype, shape.n_var)
        _append_real_field(fields, "dxTrajectoryMin", real_dtype, shape.n_var)
        _append_real_field(fields, "auxTrajectoryMax", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMin", real_dtype, shape.n_aux)
        _append_real_field(fields, "auxTrajectoryMean", real_dtype, shape.n_aux)
        _append_real_field(fields, "tUpTransition", real_dtype, n_store_events)
        _append_real_field(fields, "tDownTransition", real_dtype, n_store_events)
        _append_real_field(fields, "nMaxima", real_dtype, 3)
        _append_real_field(fields, "period", real_dtype, 3)
        _append_real_field(fields, "upDuration", real_dtype, 3)
        _append_real_field(fields, "downDuration", real_dtype, 3)
        _append_real_field(fields, "duty", real_dtype, 3)
        _append_real_field(fields, "activeDip", real_dtype, 3)
        _append_real_field(fields, "stepDt", real_dtype, 3)
        _append_real_field(fields, "fVarUpstateMean", real_dtype)
        _append_real_field(fields, "fVarDownstateMean", real_dtype)
        _append_real_field(fields, "xGlobalMax", real_dtype)
        _append_real_field(fields, "xGlobalMin", real_dtype)
        _append_real_field(fields, "dxGlobalMax", real_dtype)
        _append_real_field(fields, "dxGlobalMin", real_dtype)
        _append_real_field(fields, "xUp", real_dtype)
        _append_real_field(fields, "xDown", real_dtype)
        _append_real_field(fields, "dxUp", real_dtype)
        _append_real_field(fields, "dxDown", real_dtype)
        _append_real_field(fields, "t_start", real_dtype)
        _append_real_field(fields, "tLastEvent", real_dtype)
        _append_real_field(fields, "tThisDown", real_dtype)
        _append_real_field(fields, "tLastMax", real_dtype)
        _append_real_field(fields, "tLastMin", real_dtype)
        _append_real_field(fields, "xLastMin", real_dtype)
        _append_uint_field(fields, "thisNMaxima")
        _append_uint_field(fields, "stepcount")
        _append_uint_field(fields, "eventcount")
        _append_uint_field(fields, "inUpstate")
        return np.dtype(fields, align=True)

    raise ValueError(f"Unsupported observer metadata request: {observer_name}")


def _observer_data_struct_name(
    observer_name: str,
    shape: ProblemShape,
    precision: Precision,
    n_store_events: int,
) -> str:
    precision_name = "float" if precision is Precision.SINGLE else "double"
    return (
        f"clode_observer_data_{observer_name}_{precision_name}"
        f"_v{shape.n_var}_a{shape.n_aux}_e{n_store_events}"
    )


def _append_real_field(
    fields: list[tuple[object, ...]],
    name: str,
    real_dtype: np.dtype,
    count: int = 1,
) -> None:
    if count <= 0:
        return
    if count == 1:
        fields.append((name, real_dtype))
        return
    fields.append((name, real_dtype, (count,)))


def _append_uint_field(
    fields: list[tuple[object, ...]],
    name: str,
    count: int = 1,
) -> None:
    if count <= 0:
        return
    if count == 1:
        fields.append((name, np.uint32))
        return
    fields.append((name, np.uint32, (count,)))


def _name_at(names: list[str], index: int) -> str:
    if not names:
        raise ValueError("Observer feature naming requires at least one variable")
    clamped_index = min(max(int(index), 0), len(names) - 1)
    return names[clamped_index]