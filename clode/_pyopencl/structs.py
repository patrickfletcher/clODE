from __future__ import annotations

from dataclasses import dataclass
import warnings

import numpy as np

from ..types import ObserverParams, SolverParams
from .models import Precision
from .runtime import OpenCLRuntime, _require_pyopencl


@dataclass(frozen=True, slots=True)
class MatchedStruct:
    name: str
    dtype: np.dtype
    c_declaration: str


def _solver_params_base_dtype(precision: Precision) -> np.dtype:
    real_dtype = np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)
    return np.dtype(
        [
            ("dt", real_dtype),
            ("dtmax", real_dtype),
            ("abstol", real_dtype),
            ("reltol", real_dtype),
            ("max_steps", np.uint32),
            ("max_store", np.uint32),
            ("nout", np.uint32),
        ],
        align=True,
    )


def _observer_params_base_dtype(precision: Precision) -> np.dtype:
    real_dtype = np.dtype(np.float32 if precision is Precision.SINGLE else np.float64)
    return np.dtype(
        [
            ("eVarIx", np.uint32),
            ("fVarIx", np.uint32),
            ("maxEventCount", np.uint32),
            ("maxEventTimestamps", np.uint32),
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

    pyopencl = _require_pyopencl()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        matched_dtype, c_declaration = pyopencl.tools.match_dtype_to_c_struct(
            runtime.device, name, base_dtype
        )
    runtime.struct_cache[name] = (matched_dtype, c_declaration)
    return MatchedStruct(name, matched_dtype, c_declaration)


def match_struct_dtype(
    runtime: OpenCLRuntime, name: str, base_dtype: np.dtype
) -> MatchedStruct:
    return _match_struct(runtime, name, base_dtype)


def get_solver_params_struct(
    runtime: OpenCLRuntime, precision: Precision
) -> MatchedStruct:
    name = (
        "clode_solver_params_float"
        if precision is Precision.SINGLE
        else "clode_solver_params_double"
    )
    return _match_struct(runtime, name, _solver_params_base_dtype(precision))


def get_observer_params_struct(
    runtime: OpenCLRuntime, precision: Precision
) -> MatchedStruct:
    name = (
        "clode_observer_params_float"
        if precision is Precision.SINGLE
        else "clode_observer_params_double"
    )
    return _match_struct(runtime, name, _observer_params_base_dtype(precision))


def pack_solver_params(
    runtime: OpenCLRuntime, solver_params: SolverParams, precision: Precision
) -> np.ndarray:
    struct_spec = get_solver_params_struct(runtime, precision)
    return np.array(
        (
            solver_params.dt,
            solver_params.dtmax,
            solver_params.abstol,
            solver_params.reltol,
            solver_params.max_steps,
            solver_params.max_store,
            solver_params.nout,
        ),
        dtype=struct_spec.dtype,
    )


def pack_observer_params(
    runtime: OpenCLRuntime, observer_params: ObserverParams, precision: Precision
) -> np.ndarray:
    struct_spec = get_observer_params_struct(runtime, precision)
    return np.array(
        (
            observer_params.e_var_ix,
            observer_params.f_var_ix,
            observer_params.max_event_count,
            observer_params.max_event_timestamps,
            observer_params.min_amp,
            observer_params.min_imi,
            observer_params.nhood_radius,
            observer_params.x_up_threshold,
            observer_params.x_down_threshold,
            observer_params.dx_up_threshold,
            observer_params.dx_down_threshold,
            observer_params.eps_dx,
        ),
        dtype=struct_spec.dtype,
    )
