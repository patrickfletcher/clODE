from pathlib import Path

import numpy as np
import pyopencl
import pyopencl.tools as cl_tools
import pytest

pytest.importorskip("pyopencl")

import clode

from clode._opencl import (
    OpenCLRuntime,
    Precision,
    get_integration_settings_struct,
    get_observer_runtime_settings_struct,
    get_trajectory_output_settings_struct,
    pack_integration_settings,
    pack_observer_runtime_settings,
    pack_trajectory_output_settings,
)
from clode.observers.types import ObserverRuntimeSettings
from clode.observers.types import EventDirection
from clode.runtime import _clode_root_dir
from clode.simulation.params import IntegrationSettings, TrajectoryOutputSettings
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID


KERNEL_ROOT = Path(_clode_root_dir)


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def _build_program(
    runtime: OpenCLRuntime,
    source: str,
    *,
    extra_options: tuple[str, ...] = (),
):
    return pyopencl.Program(runtime.context, source).build(
        options=["-DCLODE_SINGLE_PRECISION", f"-I{KERNEL_ROOT}", *extra_options]
    )


def test_internal_config_structs_use_device_matched_struct_dtypes() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())

    integration_double = get_integration_settings_struct(runtime, Precision.DOUBLE)
    trajectory_output = get_trajectory_output_settings_struct(runtime)
    observer_runtime_double = get_observer_runtime_settings_struct(
        runtime, Precision.DOUBLE
    )

    assert integration_double.dtype.itemsize == 40
    assert trajectory_output.dtype.itemsize == 8
    assert observer_runtime_double.dtype.itemsize == 80
    assert "double dt;" in integration_double.c_declaration
    assert "ulong max_steps;" in integration_double.c_declaration
    assert "uint max_store;" in trajectory_output.c_declaration
    assert "uint eVarIx;" in observer_runtime_double.c_declaration


def test_pack_helpers_preserve_internal_config_values() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    integration_settings = IntegrationSettings(
        dt=0.125,
        dtmax=0.5,
        abstol=1e-6,
        reltol=1e-3,
        max_steps=(1 << 40) + 123,
    )
    trajectory_output_settings = TrajectoryOutputSettings(
        max_store=456,
        nout=7,
    )
    observer_runtime_settings = ObserverRuntimeSettings(
        e_var_ix=1,
        f_var_ix=2,
        max_event_count=30,
        event_direction=EventDirection.either,
        min_amp=0.2,
        min_imi=0.3,
        nhood_radius=0.4,
        x_up_threshold=0.5,
        x_down_threshold=0.6,
        dx_up_threshold=0.7,
        dx_down_threshold=0.8,
        eps_dx=0.9,
    )

    packed_integration = pack_integration_settings(
        runtime,
        integration_settings,
        Precision.SINGLE,
    )
    packed_trajectory_output = pack_trajectory_output_settings(
        runtime,
        trajectory_output_settings,
    )
    packed_observer_runtime = pack_observer_runtime_settings(
        runtime,
        observer_runtime_settings,
        Precision.DOUBLE,
    )

    assert float(packed_integration["dt"]) == pytest.approx(0.125)
    assert int(packed_integration["max_steps"]) == (1 << 40) + 123
    assert int(packed_trajectory_output["max_store"]) == 456
    assert int(packed_trajectory_output["nout"]) == 7
    assert int(packed_observer_runtime["eVarIx"]) == 1
    assert int(packed_observer_runtime["eventDirection"]) == int(EventDirection.either)
    assert float(packed_observer_runtime["xDownThresh"]) == pytest.approx(0.6)
    assert float(packed_observer_runtime["eps_dx"]) == pytest.approx(0.9)


def test_basic_observer_double_formula_underestimates_matched_struct_size() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    basic_observer_dtype = np.dtype(
        [
            ("xTrajectoryMax", np.float64),
            ("xTrajectoryMin", np.float64),
            ("xTrajectoryMean", np.float64),
            ("xTrajectoryIntegral", np.float64),
            ("xTrajectoryIntegralCorrection", np.float64),
            ("dxTrajectoryMax", np.float64),
            ("dxTrajectoryMin", np.float64),
            ("t_last", np.float64),
            ("t_start", np.float64),
            ("stepcount", np.uint32),
        ],
        align=True,
    )

    matched_dtype, _ = cl_tools.match_dtype_to_c_struct(
        runtime.device, "ObserverState_basic_double", basic_observer_dtype
    )

    assert matched_dtype.itemsize == 80
    assert matched_dtype.itemsize > 9 * np.dtype(np.float64).itemsize + np.dtype(np.uint32).itemsize


@pytest.mark.parametrize(
    ("observer", "observer_kwargs", "extra_options"),
    [
        (
            clode.Observer.basic_all_variables,
            {},
            ("-DUSE_OBSERVER_BASIC_ALLVAR", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=0"),
        ),
        (
            clode.Observer.local_max,
            {"observer_max_event_timestamps": 4},
            ("-DUSE_OBSERVER_LOCAL_MAX", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
        ),
        (
            clode.Observer.neighbourhood_1,
            {},
            ("-DUSE_OBSERVER_NEIGHBORHOOD_1", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=0"),
        ),
        (
            clode.Observer.neighbourhood_2,
            {"observer_max_event_timestamps": 4},
            ("-DUSE_OBSERVER_NEIGHBORHOOD_2", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
        ),
        (
            clode.Observer.threshold_crossing,
            {"observer_max_event_timestamps": 4},
            ("-DUSE_OBSERVER_THRESHOLD_CROSSING", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
        ),
        (
            clode.Observer.normalized_threshold_crossing,
            {"observer_max_event_timestamps": 4},
            ("-DUSE_OBSERVER_NORMALIZED_THRESHOLD_CROSSING", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
        ),
        (
            clode.Observer.schmitt_trigger,
            {"observer_max_event_count": 50, "observer_max_event_timestamps": 50},
            ("-DUSE_OBSERVER_SCHMITT_TRIGGER", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=50"),
        ),
        (
            clode.Observer.normalized_schmitt_trigger,
            {"observer_max_event_count": 50, "observer_max_event_timestamps": 50},
            ("-DUSE_OBSERVER_NORMALIZED_SCHMITT_TRIGGER", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=50"),
        ),
        (
            clode.Observer.threshold_2,
            {"observer_max_event_count": 50, "observer_max_event_timestamps": 50},
            ("-DUSE_OBSERVER_THRESHOLD_2", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=50"),
        ),
    ],
)
def test_observer_metadata_struct_size_matches_kernel_struct(
    observer, observer_kwargs, extra_options
) -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    simulator = clode.FeatureSimulator(
        src_file="test/van_der_pol_oscillator.cl",
        variables={"x": 0.0, "y": 1.0},
        parameters={"mu": 1.0},
        observer=observer,
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 10.0),
        **observer_kwargs,
        **_explicit_runtime_kwargs(),
    )

    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __kernel void observer_state_size(__global ulong *out) {
            out[0] = (ulong)sizeof(ObserverState);
        }
        """,
        extra_options=extra_options,
    )

    out = np.empty(1, dtype=np.uint64)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)
    program.observer_state_size(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert int(out[0]) == simulator._integrator._feature_metadata.observer_state_nbytes