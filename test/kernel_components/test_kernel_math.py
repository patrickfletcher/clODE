import math
from pathlib import Path

import numpy as np
import pytest

pyopencl = pytest.importorskip("pyopencl")

from clode._opencl import OpenCLRuntime
from clode.runtime import _clode_root_dir
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID


KERNEL_ROOT = Path(_clode_root_dir)


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def _build_program(runtime: OpenCLRuntime, source: str, *, extra_options: tuple[str, ...] = ()):
    return pyopencl.Program(runtime.context, source).build(
        options=["-DCLODE_SINGLE_PRECISION", f"-I{KERNEL_ROOT}", *extra_options]
    )


def test_compensated_sum_helper_recovers_small_term_lost_by_naive_sum() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void accumulate_values(
            __global const realtype *values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype naive = ZERO;
            realtype compensated = ZERO;
            realtype correction = ZERO;
            for (uint idx = 0; idx < n_values; ++idx) {
                naive += values[idx];
                compensatedSumAdd(&compensated, &correction, values[idx]);
            }
            out[0] = naive;
            out[1] = compensatedSumValue(compensated, correction);
        }
        """,
    )

    values = np.array([1.0e8, 1.0, -1.0e8], dtype=np.float32)
    out = np.empty(2, dtype=np.float32)
    values_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.accumulate_values(runtime.queue, (1,), None, values_buffer, out_buffer, np.uint32(len(values)))
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(0.0)
    assert out[1] == pytest.approx(1.0)


def test_basic_observer_component_kernel_reports_expected_features() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverParams TEST_PARAMS = {
            0,
            0,
            100,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_basic_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global realtype *features,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {dx_values[0]};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                dxi[0] = dx_values[idx];
                updateObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            }
            finalizeFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS, features, 0, 1);
        }
        """,
        extra_options=("-DUSE_OBSERVER_BASIC", "-DN_VAR=1", "-DN_AUX=0"),
    )

    times = np.array([0.0, 1.0, 3.0], dtype=np.float32)
    x_values = np.array([0.0, 10.0, 4.0], dtype=np.float32)
    dx_values = np.array([1.0, 2.0, -1.0], dtype=np.float32)
    features = np.empty(6, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=times,
    )
    x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=x_values,
    )
    dx_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=dx_values,
    )
    feature_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.WRITE_ONLY,
        features.nbytes,
    )

    program.run_basic_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        dx_buffer,
        feature_buffer,
        np.uint32(len(times)),
    )
    pyopencl.enqueue_copy(runtime.queue, features, feature_buffer).wait()

    assert features[0] == pytest.approx(10.0)
    assert features[1] == pytest.approx(4.0)
    assert features[2] == pytest.approx(6.0)
    assert features[3] == pytest.approx(2.0)
    assert features[4] == pytest.approx(-1.0)
    assert features[5] == pytest.approx(2.0)


def test_basicall_observer_component_kernel_reports_expected_feature_layout() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverParams TEST_PARAMS = {
            0,
            0,
            100,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_basicall_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global const realtype *aux_values,
            __global realtype *features,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[2] = {x_values[0], x_values[1]};
            realtype dxi[2] = {dx_values[0], dx_values[1]};
            realtype auxi[1] = {aux_values[0]};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                for (uint j = 0; j < 2; ++j) {
                    xi[j] = x_values[idx * 2 + j];
                    dxi[j] = dx_values[idx * 2 + j];
                }
                auxi[0] = aux_values[idx];
                updateObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            }
            finalizeFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS, features, 0, 1);
        }
        """,
        extra_options=("-DUSE_OBSERVER_BASIC_ALLVAR", "-DN_VAR=2", "-DN_AUX=1"),
    )

    times = np.array([0.0, 1.0, 3.0], dtype=np.float32)
    x_values = np.array([0.0, 1.0, 10.0, -5.0, 4.0, 8.0], dtype=np.float32)
    dx_values = np.array([1.0, 0.0, 2.0, -1.0, -1.0, 4.0], dtype=np.float32)
    aux_values = np.array([7.0, 9.0, 6.0], dtype=np.float32)
    features = np.empty(14, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=times,
    )
    x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=x_values,
    )
    dx_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=dx_values,
    )
    aux_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=aux_values,
    )
    feature_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.WRITE_ONLY,
        features.nbytes,
    )

    program.run_basicall_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        dx_buffer,
        aux_buffer,
        feature_buffer,
        np.uint32(len(times)),
    )
    pyopencl.enqueue_copy(runtime.queue, features, feature_buffer).wait()

    expected = np.array(
        [10.0, 4.0, 6.0, 2.0, -1.0, 8.0, -5.0, 11.0 / 3.0, 4.0, -1.0, 9.0, 6.0, 7.0, 2.0],
        dtype=np.float32,
    )
    np.testing.assert_allclose(features, expected, rtol=1e-6, atol=1e-6)


def test_threshold_2_observer_kernel_uses_inverse_linear_threshold_timestamps() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverParams TEST_PARAMS = {
            0,
            0,
            4,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_threshold_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {dx_values[0]};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            observer_state.xGlobalMax = ONE;
            observer_state.xGlobalMin = -ONE;
            observer_state.xUp = RCONST(0.5);
            observer_state.xDown = -RCONST(0.5);
            observer_state.dxUp = ZERO;
            observer_state.dxDown = ZERO;
            observer_state.inUpstate = 0;

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                dxi[0] = dx_values[idx];
                updateObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tUpTransition[0];
            out[1] = observer_state.tDownTransition[0];
            out[2] = observer_state.tUpTransition[1];
            out[3] = observer_state.tDownTransition[1];
        }
        """,
        extra_options=("-DUSE_OBSERVER_THRESHOLD_2", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    dt = 0.2
    sample_times = np.arange(0.0, 4.0 * math.pi + dt, dt, dtype=np.float64)
    x_values = np.sin(sample_times).astype(np.float32)
    dx_values = np.cos(sample_times).astype(np.float32)
    out = np.empty(4, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=sample_times.astype(np.float32),
    )
    x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=x_values,
    )
    dx_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=dx_values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_threshold_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        dx_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(
        out,
        np.array([math.pi / 6.0, 7.0 * math.pi / 6.0, 13.0 * math.pi / 6.0, 19.0 * math.pi / 6.0], dtype=np.float32),
        atol=2e-2,
        rtol=0.0,
    )


def test_threshold_2_observer_kernel_ignores_dx_when_dx_thresholds_zero() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverParams TEST_PARAMS = {
            0,
            0,
            4,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_threshold_observer_ignore_dx(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {dx_values[0]};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            observer_state.xGlobalMax = ONE;
            observer_state.xGlobalMin = -ONE;
            observer_state.xUp = RCONST(0.5);
            observer_state.xDown = -RCONST(0.5);
            observer_state.dxUp = RCONST(0.25);
            observer_state.dxDown = -RCONST(0.25);
            observer_state.inUpstate = 0;

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                dxi[0] = dx_values[idx];
                updateObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tUpTransition[0];
            out[1] = observer_state.tDownTransition[0];
            out[2] = observer_state.eventcount;
        }
        """,
        extra_options=("-DUSE_OBSERVER_THRESHOLD_2", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float32)
    x_values = np.array([-1.0, 0.0, 1.0, 0.0, -1.0], dtype=np.float32)
    dx_values = np.array([-1.0, -1.0, -1.0, 1.0, 1.0], dtype=np.float32)
    out = np.empty(3, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=sample_times,
    )
    x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=x_values,
    )
    dx_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=dx_values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_threshold_observer_ignore_dx(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        dx_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(out[:2], np.array([1.5, 3.5], dtype=np.float32), atol=1e-6, rtol=0.0)
    assert out[2] == 1.0


@pytest.mark.parametrize(
    ("name", "rising", "x0", "x1", "dx0", "dx1", "x_threshold", "dx_threshold", "expected"),
    [
        ("rising_x_first", True, 0.0, 1.0, 0.0, 1.0, 0.25, 0.75, 0.75),
        ("rising_dx_first", True, 0.0, 1.0, 0.0, 1.0, 0.75, 0.25, 0.75),
        ("falling_x_first", False, 0.0, -1.0, 0.0, -1.0, -0.25, -0.75, 0.75),
        ("falling_dx_first", False, 0.0, -1.0, 0.0, -1.0, -0.75, -0.25, 0.75),
    ],
)
def test_threshold_transition_time_returns_later_active_boundary(
    name: str,
    rising: bool,
    x0: float,
    x1: float,
    dx0: float,
    dx1: float,
    x_threshold: float,
    dx_threshold: float,
    expected: float,
) -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __kernel void evaluate_transition_time(
            const uint use_dx,
            const uint rising,
            const realtype t0,
            const realtype t1,
            const realtype x0,
            const realtype x1,
            const realtype dx0,
            const realtype dx1,
            const realtype x_threshold,
            const realtype dx_threshold,
            __global realtype *out
        ) {
            out[0] = thresholdTransitionTime(
                use_dx != 0,
                rising != 0,
                t0,
                t1,
                x0,
                x1,
                dx0,
                dx1,
                x_threshold,
                dx_threshold
            );
        }
        """,
        extra_options=("-DUSE_OBSERVER_THRESHOLD_2", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    out = np.empty(1, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.evaluate_transition_time(
        runtime.queue,
        (1,),
        None,
        np.uint32(1),
        np.uint32(1 if rising else 0),
        np.float32(0.0),
        np.float32(1.0),
        np.float32(x0),
        np.float32(x1),
        np.float32(dx0),
        np.float32(dx1),
        np.float32(x_threshold),
        np.float32(dx_threshold),
        out_buffer,
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(expected, abs=1e-6), name


@pytest.mark.parametrize(
    ("name", "x_values", "dx_values", "expected_up", "expected_down"),
    [
        (
            "x_first_dx_later",
            np.array([-1.0, 0.6, 1.0, -0.6, -1.0], dtype=np.float32),
            np.array([-1.0, 0.0, 0.5, 0.0, -0.5], dtype=np.float32),
            1.5,
            3.5,
        ),
        (
            "dx_first_x_later",
            np.array([-1.0, 0.0, 0.6, 0.0, -0.6], dtype=np.float32),
            np.array([-1.0, 0.5, 1.0, -0.5, -1.0], dtype=np.float32),
            1.0 + 0.5 / 0.6,
            3.0 + 0.5 / 0.6,
        ),
    ],
)
def test_threshold_2_state_machine_waits_for_second_active_gate(
    name: str,
    x_values: np.ndarray,
    dx_values: np.ndarray,
    expected_up: float,
    expected_down: float,
) -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverParams TEST_PARAMS = {
            0,
            0,
            4,
            ZERO,
            ZERO,
            ZERO,
            RCONST(0.75),
            RCONST(0.25),
            RCONST(0.25),
            RCONST(0.25)
        };

        __kernel void run_threshold_observer_waits_for_second_gate(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {dx_values[0]};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            observer_state.xGlobalMax = ONE;
            observer_state.xGlobalMin = -ONE;
            observer_state.dxGlobalMax = ONE;
            observer_state.dxGlobalMin = -ONE;
            initializeEventDetector(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                dxi[0] = dx_values[idx];
                updateObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tUpTransition[0];
            out[1] = observer_state.tDownTransition[0];
            out[2] = observer_state.eventcount;
        }
        """,
        extra_options=("-DUSE_OBSERVER_THRESHOLD_2", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.array([0.0, 1.0, 2.0, 3.0, 4.0], dtype=np.float32)
    out = np.empty(3, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=sample_times,
    )
    x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=x_values,
    )
    dx_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=dx_values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_threshold_observer_waits_for_second_gate(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        dx_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(out[:2], np.array([expected_up, expected_down], dtype=np.float32), atol=1e-6, rtol=0.0)
    assert out[2] == pytest.approx(1.0), name


def test_local_max_observer_kernel_uses_three_sample_extremum_helpers() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverParams TEST_PARAMS = {
            0,
            0,
            4,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_localmax_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {dx_values[0]};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                dxi[0] = dx_values[idx];
                updateObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tMaxList[0];
            out[1] = observer_state.tMinList[0];
            out[2] = observer_state.tMaxList[1];
            out[3] = observer_state.tMinList[1];
        }
        """,
        extra_options=("-DUSE_OBSERVER_LOCAL_MAX", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    dt = 0.2
    sample_times = np.arange(0.0, 4.0 * math.pi + dt, dt, dtype=np.float64)
    x_values = np.sin(sample_times).astype(np.float32)
    dx_values = np.cos(sample_times).astype(np.float32)
    out = np.empty(4, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=sample_times.astype(np.float32),
    )
    x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=x_values,
    )
    dx_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=dx_values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_localmax_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        dx_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(
        out,
        np.array([math.pi / 2.0, 3.0 * math.pi / 2.0, 5.0 * math.pi / 2.0, 7.0 * math.pi / 2.0], dtype=np.float32),
        atol=2e-2,
        rtol=0.0,
    )


def test_relative_elapsed_mean_prototype_survives_large_origin_bias() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void compare_elapsed_strategies(
            __global const realtype *values,
            __global realtype *out,
            const uint n_values
        ) {
            const realtype t0 = RCONST(10000.0);
            const realtype dt = RCONST(0.01);
            realtype ti = t0;
            realtype t_start = t0;
            realtype integral = ZERO;
            realtype integral_correction = ZERO;
            realtype elapsed = ZERO;
            realtype elapsed_correction = ZERO;

            for (uint idx = 0; idx < n_values; ++idx) {
                ti += dt;
                compensatedIntegrateConstant(&integral, &integral_correction, dt, values[idx]);
                compensatedTimeAdd(&elapsed, &elapsed_correction, dt);
            }

            out[0] = meanFromCompensatedIntegral(integral, integral_correction, ti - t_start);
            out[1] = meanFromCompensatedIntegral(
                integral,
                integral_correction,
                compensatedTimeValue(elapsed, elapsed_correction)
            );
            out[2] = ti - t_start;
            out[3] = compensatedTimeValue(elapsed, elapsed_correction);
        }
        """,
    )

    values = np.ones(10_000, dtype=np.float32)
    values[5_000:] += np.float32(1.0e-4)
    expected = float(values.mean(dtype=np.float64))
    out = np.empty(4, dtype=np.float32)

    values_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.compare_elapsed_strategies(
        runtime.queue,
        (1,),
        None,
        values_buffer,
        out_buffer,
        np.uint32(len(values)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] > np.float32(1.02)
    assert out[1] == pytest.approx(expected, abs=2e-6)
    assert out[2] == pytest.approx(np.float32(97.65625))
    assert out[3] == pytest.approx(np.float32(100.0))


def test_time_prototypes_reduce_large_origin_failure_and_zero_origin_drift() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void compare_time_prototypes(__global realtype *out) {
            const realtype dt = RCONST(0.01);

            realtype large_direct = RCONST(1000000.0);
            realtype large_kahan = RCONST(1000000.0);
            realtype large_kahan_correction = ZERO;
            for (uint step = 0; step < 10000; ++step) {
                large_direct += dt;
                compensatedTimeAdd(&large_kahan, &large_kahan_correction, dt);
            }
            out[0] = large_direct;
            out[1] = large_kahan;
            out[2] = large_kahan_correction;
            out[3] = fixedStepTimeFromCounter(RCONST(1000000.0), (ulong)10000, dt);

            realtype zero_direct = ZERO;
            realtype zero_kahan = ZERO;
            realtype zero_kahan_correction = ZERO;
            for (uint step = 0; step < 200000; ++step) {
                zero_direct += dt;
                compensatedTimeAdd(&zero_kahan, &zero_kahan_correction, dt);
            }
            out[4] = zero_direct;
            out[5] = zero_kahan;
            out[6] = zero_kahan_correction;
            out[7] = fixedStepTimeFromCounter(ZERO, (ulong)200000, dt);

            realtype step3 = fixedStepTimeFromCounter(RCONST(1000000.0), (ulong)3, dt);
            realtype step4 = fixedStepTimeFromCounter(RCONST(1000000.0), (ulong)4, dt);
            out[8] = step3 + RCONST(0.5) * (step4 - step3);
            out[9] = fixedStepTimeFromRealIndex(RCONST(1000000.0), RCONST(3.5), dt);
        }
        """,
    )

    out = np.empty(10, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.compare_time_prototypes(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(np.float32(1.0e6))
    assert out[1] == pytest.approx(np.float32(1000100.0))
    assert out[3] == pytest.approx(np.float32(1000100.0))
    assert abs(float(out[4]) - 2000.0) > 1.0
    assert out[5] == pytest.approx(np.float32(2000.0))
    assert out[7] == pytest.approx(np.float32(2000.0))
    assert out[8] == pytest.approx(np.float32(1000000.0))
    assert out[9] == pytest.approx(np.float32(1000000.0625))
    assert out[9] > out[8]


def test_threshold_crossing_interpolation_helpers_improve_timestamp_accuracy() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void interpolate_threshold_crossing(
            __global const realtype *samples,
            __global realtype *out
        ) {
            realtype t0 = samples[0];
            realtype t1 = samples[1];
            realtype x0 = samples[2];
            realtype x1 = samples[3];
            realtype dx0 = samples[4];
            realtype dx1 = samples[5];
            realtype threshold = samples[6];

            out[0] = t1;
            out[1] = linearInterpTimeOfValue(t0, t1, x0, x1, threshold);
            out[2] = cubicHermiteInterpTimeOfValue(t0, t1, x0, x1, dx0, dx1, threshold);
        }
        """,
    )

    dt = 0.2
    offset = 0.03
    threshold = 0.5
    sample_times = np.arange(0.0, 8.0 * math.pi, dt, dtype=np.float64) + offset
    x = np.sin(sample_times)
    dx = np.cos(sample_times)
    crossing_index = next(
        index
        for index in range(1, len(sample_times))
        if x[index - 1] <= threshold < x[index]
    )
    true_time = math.asin(threshold)
    samples = np.array(
        [
            sample_times[crossing_index - 1],
            sample_times[crossing_index],
            x[crossing_index - 1],
            x[crossing_index],
            dx[crossing_index - 1],
            dx[crossing_index],
            threshold,
        ],
        dtype=np.float32,
    )
    out = np.empty(3, dtype=np.float32)

    samples_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=samples,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.interpolate_threshold_crossing(runtime.queue, (1,), None, samples_buffer, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    sample_error = abs(float(out[0]) - true_time)
    linear_error = abs(float(out[1]) - true_time)
    hermite_error = abs(float(out[2]) - true_time)

    assert linear_error < sample_error / 10.0
    assert hermite_error < linear_error / 100.0


def test_three_sample_extremum_helpers_improve_over_sample_pick() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void recover_local_extrema(
            __global const realtype *max_t,
            __global const realtype *max_x,
            __global const realtype *min_t,
            __global const realtype *min_x,
            __global realtype *out
        ) {
            realtype max_time_buffer[3] = {max_t[0], max_t[1], max_t[2]};
            realtype max_value_buffer[3] = {max_x[0], max_x[1], max_x[2]};
            realtype min_time_buffer[3] = {min_t[0], min_t[1], min_t[2]};
            realtype min_value_buffer[3] = {min_x[0], min_x[1], min_x[2]};

            int max_index = array_argmax(max_value_buffer, 3);
            int min_index = array_argmin(min_value_buffer, 3);

            out[0] = max_time_buffer[max_index];
            out[1] = max_value_buffer[max_index];
            localMaximumFromThreeSamples(max_time_buffer, max_value_buffer, &out[2], &out[3]);

            out[4] = min_time_buffer[min_index];
            out[5] = min_value_buffer[min_index];
            localMinimumFromThreeSamples(min_time_buffer, min_value_buffer, &out[6], &out[7]);
        }
        """,
    )

    dt = 0.2
    offset = 0.03
    sample_times = np.arange(0.0, 12.0 * math.pi, dt, dtype=np.float64) + offset
    x = np.sin(sample_times)
    dx = np.cos(sample_times)

    max_index = next(index for index in range(2, len(sample_times)) if dx[index - 1] > 0.0 and dx[index] < 0.0)
    min_index = next(index for index in range(2, len(sample_times)) if dx[index - 1] < 0.0 and dx[index] > 0.0)

    max_t = sample_times[max_index - 2 : max_index + 1].astype(np.float32)
    max_x = x[max_index - 2 : max_index + 1].astype(np.float32)
    min_t = sample_times[min_index - 2 : min_index + 1].astype(np.float32)
    min_x = x[min_index - 2 : min_index + 1].astype(np.float32)
    out = np.empty(8, dtype=np.float32)

    max_t_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=max_t,
    )
    max_x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=max_x,
    )
    min_t_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=min_t,
    )
    min_x_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=min_x,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.recover_local_extrema(
        runtime.queue,
        (1,),
        None,
        max_t_buffer,
        max_x_buffer,
        min_t_buffer,
        min_x_buffer,
        out_buffer,
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    max_sample_error = abs(float(out[0]) - math.pi / 2.0)
    max_helper_error = abs(float(out[2]) - math.pi / 2.0)
    min_sample_error = abs(float(out[4]) - 3.0 * math.pi / 2.0)
    min_helper_error = abs(float(out[6]) - 3.0 * math.pi / 2.0)

    assert max_helper_error < max_sample_error / 50.0
    assert min_helper_error < min_sample_error / 50.0
    assert abs(float(out[3]) - 1.0) < abs(float(out[1]) - 1.0)
    assert abs(float(out[7]) + 1.0) < abs(float(out[5]) + 1.0)