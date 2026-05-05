import numpy as np
import pyopencl.tools as cl_tools
import pytest

pytest.importorskip("pyopencl")

from clode._pyopencl import (
    OpenCLRuntime,
    Precision,
    get_observer_params_struct,
    get_solver_params_struct,
    pack_observer_params,
    pack_solver_params,
)
from clode.cpp.clode_cpp_wrapper import ObserverParams, SolverParams
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def test_solver_and_observer_params_use_device_matched_struct_dtypes() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())

    solver_double = get_solver_params_struct(runtime, Precision.DOUBLE)
    observer_double = get_observer_params_struct(runtime, Precision.DOUBLE)

    assert solver_double.dtype.itemsize == 48
    assert observer_double.dtype.itemsize == 80
    assert "double dt;" in solver_double.c_declaration
    assert "uint eVarIx;" in observer_double.c_declaration


def test_pack_helpers_preserve_solver_and_observer_param_values() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    solver_params = SolverParams(
        dt=0.125,
        dtmax=0.5,
        abstol=1e-6,
        reltol=1e-3,
        max_steps=123,
        max_store=456,
        nout=7,
    )
    observer_params = ObserverParams(
        e_var_ix=1,
        f_var_ix=2,
        max_event_count=30,
        max_event_timestamps=4,
        min_amp=0.2,
        min_imi=0.3,
        nhood_radius=0.4,
        x_up_threshold=0.5,
        x_down_threshold=0.6,
        dx_up_threshold=0.7,
        dx_down_threshold=0.8,
        eps_dx=0.9,
    )

    packed_solver = pack_solver_params(runtime, solver_params, Precision.SINGLE)
    packed_observer = pack_observer_params(runtime, observer_params, Precision.DOUBLE)

    assert float(packed_solver["dt"]) == pytest.approx(0.125)
    assert int(packed_solver["max_steps"]) == 123
    assert int(packed_observer["eVarIx"]) == 1
    assert float(packed_observer["xDownThresh"]) == pytest.approx(0.6)
    assert float(packed_observer["eps_dx"]) == pytest.approx(0.9)


def test_basic_observer_double_formula_underestimates_matched_struct_size() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    basic_observer_dtype = np.dtype(
        [
            ("xTrajectoryMax", np.float64),
            ("xTrajectoryMin", np.float64),
            ("xTrajectoryMean", np.float64),
            ("dxTrajectoryMax", np.float64),
            ("dxTrajectoryMin", np.float64),
            ("t_last", np.float64),
            ("t_start", np.float64),
            ("stepcount", np.uint32),
        ],
        align=True,
    )

    matched_dtype, _ = cl_tools.match_dtype_to_c_struct(
        runtime.device, "ObserverData_basic_double", basic_observer_dtype
    )

    assert matched_dtype.itemsize == 64
    assert matched_dtype.itemsize > 7 * np.dtype(np.float64).itemsize + np.dtype(np.uint32).itemsize