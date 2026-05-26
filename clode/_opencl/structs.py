from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np

from ..observers.types import ObserverParams, ObserverRuntimeSettings
from ..simulation.params import (
    SolverParams,
    IntegrationSettings,
    TrajectoryOutputSettings,
)
from .models import Precision
from .runtime import OpenCLRuntime, _require_opencl_binding


@dataclass(frozen=True, slots=True)
class MatchedStruct:
    name: str
    dtype: np.dtype
    c_declaration: str


def _integration_settings_base_dtype(precision: Precision) -> np.dtype:
    real_dtype = np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)
    return np.dtype(
        [
            ("dt", real_dtype),
            ("dtmax", real_dtype),
            ("abstol", real_dtype),
            ("reltol", real_dtype),
            ("max_steps", np.uint64),
        ],
        align=True,
    )


def _trajectory_output_settings_base_dtype() -> np.dtype:
    return np.dtype(
        [
            ("max_store", np.uint32),
            ("nout", np.uint32),
        ],
        align=True,
    )


def _observer_runtime_settings_base_dtype(precision: Precision) -> np.dtype:
    real_dtype = np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)
    return np.dtype(
        [
            ("eVarIx", np.uint32),
            ("fVarIx", np.uint32),
            ("maxEventCount", np.uint32),
            ("minXamp", real_dtype),
            ("minIMI", real_dtype),
            ("nHoodRadius", real_dtype),
            ("xUpThresh", real_dtype),
            ("xDownThresh", real_dtype),
            ("dxUpThresh", real_dtype),
            ("dxDownThresh", real_dtype),
            ("eps_dx", real_dtype),
        ],
        align=True,
    )


def _match_struct(
    runtime: OpenCLRuntime, name: str, base_dtype: np.dtype
) -> MatchedStruct:
    cached = runtime.struct_cache.get(name)
    if cached is not None:
        return MatchedStruct(name, cached[0], cached[1])

    opencl_binding = _require_opencl_binding()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        matched_dtype, c_declaration = opencl_binding.tools.match_dtype_to_c_struct(
            runtime.device, name, base_dtype
        )
    runtime.struct_cache[name] = (matched_dtype, c_declaration)
    return MatchedStruct(name, matched_dtype, c_declaration)


def match_struct_dtype(
    runtime: OpenCLRuntime, name: str, base_dtype: np.dtype
) -> MatchedStruct:
    return _match_struct(runtime, name, base_dtype)


def get_integration_settings_struct(
    runtime: OpenCLRuntime, precision: Precision
) -> MatchedStruct:
    name = (
        "clode_integration_settings_float"
        if precision is Precision.SINGLE
        else "clode_integration_settings_double"
    )
    return _match_struct(runtime, name, _integration_settings_base_dtype(precision))


def get_trajectory_output_settings_struct(runtime: OpenCLRuntime) -> MatchedStruct:
    return _match_struct(
        runtime,
        "clode_trajectory_output_settings",
        _trajectory_output_settings_base_dtype(),
    )


def get_observer_runtime_settings_struct(
    runtime: OpenCLRuntime, precision: Precision
) -> MatchedStruct:
    name = (
        "clode_observer_runtime_settings_float"
        if precision is Precision.SINGLE
        else "clode_observer_runtime_settings_double"
    )
    return _match_struct(runtime, name, _observer_runtime_settings_base_dtype(precision))


def pack_integration_settings(
    runtime: OpenCLRuntime,
    integration_settings: IntegrationSettings,
    precision: Precision,
) -> np.ndarray:
    struct_spec = get_integration_settings_struct(runtime, precision)
    return np.array(
        (
            integration_settings.dt,
            integration_settings.dtmax,
            integration_settings.abstol,
            integration_settings.reltol,
            integration_settings.max_steps,
        ),
        dtype=struct_spec.dtype,
    )


def pack_trajectory_output_settings(
    runtime: OpenCLRuntime,
    output_settings: TrajectoryOutputSettings,
) -> np.ndarray:
    struct_spec = get_trajectory_output_settings_struct(runtime)
    return np.array(
        (
            output_settings.max_store,
            output_settings.nout,
        ),
        dtype=struct_spec.dtype,
    )


def pack_observer_runtime_settings(
    runtime: OpenCLRuntime,
    observer_runtime_settings: ObserverRuntimeSettings,
    precision: Precision,
) -> np.ndarray:
    struct_spec = get_observer_runtime_settings_struct(runtime, precision)
    return np.array(
        (
            observer_runtime_settings.e_var_ix,
            observer_runtime_settings.f_var_ix,
            observer_runtime_settings.max_event_count,
            observer_runtime_settings.min_amp,
            observer_runtime_settings.min_imi,
            observer_runtime_settings.nhood_radius,
            observer_runtime_settings.x_up_threshold,
            observer_runtime_settings.x_down_threshold,
            observer_runtime_settings.dx_up_threshold,
            observer_runtime_settings.dx_down_threshold,
            observer_runtime_settings.eps_dx,
        ),
        dtype=struct_spec.dtype,
    )


def pack_solver_params(
    runtime: OpenCLRuntime, solver_params: SolverParams, precision: Precision
) -> np.ndarray:
    return pack_integration_settings(
        runtime,
        solver_params.integration_settings,
        precision,
    )


def pack_observer_params(
    runtime: OpenCLRuntime, observer_params: ObserverParams, precision: Precision
) -> np.ndarray:
    return pack_observer_runtime_settings(
        runtime,
        observer_params.runtime_settings,
        precision,
    )


def get_solver_params_struct(
    runtime: OpenCLRuntime, precision: Precision
) -> MatchedStruct:
    return get_integration_settings_struct(runtime, precision)


def get_observer_params_struct(
    runtime: OpenCLRuntime, precision: Precision
) -> MatchedStruct:
    return get_observer_runtime_settings_struct(runtime, precision)
