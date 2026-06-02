import math
from pathlib import Path

import numpy as np
import pytest

pyopencl = pytest.importorskip("pyopencl")

from clode._opencl import OpenCLRuntime
from clode.observers._definitions import resolve_observer_spec
from clode.observers.types import EventOutputSettings
from clode.observers.types import EventDirection
from clode.observers.types import ObserverRuntimeSettings
from clode.problem._core import ProblemInfo
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


def _summary_observer_source_preamble(
    observer_name: str,
    *,
    n_var: int,
    n_aux: int,
) -> str:
    problem_info = ProblemInfo(
        "component_kernel.cl",
        [f"x{index}" for index in range(n_var)],
        ["p0"],
        [f"aux{index}" for index in range(n_aux)],
        0,
    )
    resolved_spec = resolve_observer_spec(
        problem_info,
        observer_name,
        ObserverRuntimeSettings(),
        real_dtype=np.dtype(np.float32),
        event_output_settings=EventOutputSettings(max_event_timestamps=0),
    )
    return resolved_spec.source_preamble


def _reference_compensated_time_add(
    value: np.float32,
    correction: np.float32,
    dt: np.float32,
) -> tuple[np.float32, np.float32]:
    y = np.float32(dt - correction)
    total = np.float32(value + y)
    correction = np.float32((total - value) - y)
    value = total
    return value, correction


def _reference_compensated_time_value(value: np.float32, correction: np.float32) -> np.float32:
    return np.float32(value - correction)


def _reference_compensated_time_from_origin(
    origin: np.float32,
    elapsed: np.float32,
    correction: np.float32,
) -> np.float32:
    return np.float32(origin + _reference_compensated_time_value(elapsed, correction))


def _reference_compensated_time_from_origin_after_step(
    origin: np.float32,
    elapsed: np.float32,
    correction: np.float32,
    dt: np.float32,
) -> np.float32:
    elapsed, correction = _reference_compensated_time_add(elapsed, correction, dt)
    return _reference_compensated_time_from_origin(origin, elapsed, correction)


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


def test_named_address_space_helpers_cover_private_and_global_array_patterns() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"realtype.cl\"

        typedef struct {
            realtype values[4];
        } HelperState;

        static inline realtype sum_private(__private realtype values[], const int n_values) {
            realtype total = ZERO;
            for (int idx = 0; idx < n_values; ++idx)
                total += values[idx];
            return total;
        }

        static inline realtype sum_global(__global const realtype *values, const int n_values) {
            realtype total = ZERO;
            for (int idx = 0; idx < n_values; ++idx)
                total += values[idx];
            return total;
        }

        __kernel void probe_named_address_space_helpers(
            __global const realtype *global_values,
            __global realtype *out
        ) {
            realtype stack_values[4] = {
                RCONST(1.0),
                RCONST(2.0),
                RCONST(3.0),
                RCONST(4.0)
            };
            HelperState state = {{
                RCONST(5.0),
                RCONST(6.0),
                RCONST(7.0),
                RCONST(8.0)
            }};

            out[0] = sum_private(stack_values, 4);
            out[1] = sum_private(state.values, 4);
            out[2] = sum_global(global_values, 4);
        }
        """,
    )

    global_values = np.array([9.0, 10.0, 11.0, 12.0], dtype=np.float32)
    out = np.empty(3, dtype=np.float32)
    global_values_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=global_values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.probe_named_address_space_helpers(
        runtime.queue,
        (1,),
        None,
        global_values_buffer,
        out_buffer,
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(
        out,
        np.array([10.0, 26.0, 42.0], dtype=np.float32),
        atol=0.0,
        rtol=0.0,
    )


def test_basic_observer_component_kernel_reports_expected_features() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    source_preamble = _summary_observer_source_preamble("summary", n_var=1, n_aux=0)
    program = _build_program(
        runtime,
        source_preamble
        + """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
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
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
            }
            finalizeFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS, times[n_values - 1] - times[0], features, 0, 1);
        }
        """,
        extra_options=("-DUSE_OBSERVER_SUMMARY", "-DN_VAR=1", "-DN_AUX=0"),
    )

    times = np.array([0.0, 1.0, 3.0], dtype=np.float32)
    x_values = np.array([0.0, 10.0, 4.0], dtype=np.float32)
    dx_values = np.array([1.0, 2.0, -1.0], dtype=np.float32)
    features = np.empty(5, dtype=np.float32)

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


def test_basic_observer_constant_signal_preserves_exact_unit_mean() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    source_preamble = _summary_observer_source_preamble("summary", n_var=1, n_aux=0)
    program = _build_program(
        runtime,
        source_preamble
        + """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
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

        __kernel void run_basic_constant_observer(
            __global const realtype *times,
            __global realtype *features,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {ONE};
            realtype dxi[1] = {ZERO};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
            }
            finalizeFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS, times[n_values - 1] - times[0], features, 0, 1);
        }
        """,
        extra_options=("-DUSE_OBSERVER_SUMMARY", "-DN_VAR=1", "-DN_AUX=0"),
    )

    times = np.arange(0.0, 400.0 * math.pi + 0.1, 0.1, dtype=np.float32)
    features = np.empty(5, dtype=np.float32)

    time_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=times,
    )
    feature_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.WRITE_ONLY,
        features.nbytes,
    )

    program.run_basic_constant_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        feature_buffer,
        np.uint32(len(times)),
    )
    pyopencl.enqueue_copy(runtime.queue, features, feature_buffer).wait()

    assert features[0] == pytest.approx(1.0)
    assert features[1] == pytest.approx(1.0)
    assert features[2] == 1.0
    assert features[3] == pytest.approx(0.0)
    assert features[4] == pytest.approx(0.0)


def test_summary_observer_component_kernel_reports_expected_feature_layout() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    source_preamble = _summary_observer_source_preamble("summary", n_var=2, n_aux=1)
    program = _build_program(
        runtime,
        source_preamble
        + """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
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

        __kernel void run_summary_observer(
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
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
            }
            finalizeFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS, times[n_values - 1] - times[0], features, 0, 1);
        }
        """,
        extra_options=("-DUSE_OBSERVER_SUMMARY", "-DN_VAR=2", "-DN_AUX=1"),
    )

    times = np.array([0.0, 1.0, 3.0], dtype=np.float32)
    x_values = np.array([0.0, 1.0, 10.0, -5.0, 4.0, 8.0], dtype=np.float32)
    dx_values = np.array([1.0, 0.0, 2.0, -1.0, -1.0, 4.0], dtype=np.float32)
    aux_values = np.array([7.0, 9.0, 6.0], dtype=np.float32)
    features = np.empty(13, dtype=np.float32)

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

    program.run_summary_observer(
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
        [10.0, 4.0, 6.0, 2.0, -1.0, 8.0, -5.0, 11.0 / 3.0, 4.0, -1.0, 9.0, 6.0, 7.0],
        dtype=np.float32,
    )
    np.testing.assert_allclose(features, expected, rtol=1e-6, atol=1e-6)



@pytest.mark.parametrize(
    ("build_define", "x_up_threshold", "x_down_threshold", "needs_warmup"),
    [
        ("USE_OBSERVER_SCHMITT_TRIGGER", 0.5, -0.5, False),
        ("USE_OBSERVER_NORMALIZED_SCHMITT_TRIGGER", 0.75, 0.25, True),
    ],
)
def test_schmitt_observer_kernels_store_up_and_down_transitions_in_event_features(
    build_define: str,
    x_up_threshold: float,
    x_down_threshold: float,
    needs_warmup: bool,
) -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        (
            """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
            0,
            0,
            4,
            ZERO,
            RCONST(0.5),
            ZERO,
            ZERO,
            __X_UP__,
            __X_DOWN__,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_schmitt_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {ZERO};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            if (__NEEDS_WARMUP__) {
                observer_state.xGlobalMax = ONE;
                observer_state.xGlobalMin = -ONE;
            }
            initializeEventDetector(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tUpTransition[0];
            out[1] = observer_state.xUpTransition[0];
            out[2] = observer_state.tDownTransition[0];
            out[3] = observer_state.xDownTransition[0];
            out[4] = observer_state.tUpTransition[1];
            out[5] = observer_state.xUpTransition[1];
            out[6] = observer_state.tDownTransition[1];
            out[7] = observer_state.xDownTransition[1];
            out[8] = (realtype)observer_state.eventcount;
        }
        """
        )
        .replace("__X_UP__", f"RCONST({x_up_threshold})")
        .replace("__X_DOWN__", f"RCONST({x_down_threshold})")
        .replace("__NEEDS_WARMUP__", "1" if needs_warmup else "0"),
        extra_options=(f"-D{build_define}", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0], dtype=np.float32)
    x_values = np.array([-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0], dtype=np.float32)
    out = np.empty(9, dtype=np.float32)

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
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_schmitt_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    expected_up_value = 0.5 if needs_warmup else x_up_threshold
    expected_down_value = -0.5 if needs_warmup else x_down_threshold
    np.testing.assert_allclose(
        out[:8],
        np.array([1.5, expected_up_value, 3.5, expected_down_value, 5.5, expected_up_value, 7.5, expected_down_value], dtype=np.float32),
        atol=1e-6,
        rtol=0.0,
    )
    assert out[8] == pytest.approx(2.0)


@pytest.mark.parametrize(
    ("event_direction", "expected_times"),
    [
        (EventDirection.rising, np.array([1.5, 5.5], dtype=np.float32)),
        (EventDirection.falling, np.array([2.5], dtype=np.float32)),
        (EventDirection.either, np.array([1.5, 2.5, 5.5], dtype=np.float32)),
    ],
)
def test_threshold_crossing_observer_kernel_tracks_absolute_crossings(
    event_direction: EventDirection,
    expected_times: np.ndarray,
) -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        (
            """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
            0,
            0,
            8,
            __EVENT_DIRECTION__,
            ZERO,
            ZERO,
            ZERO,
            RCONST(0.5),
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_threshold1_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {ZERO};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            for (int idx = 0; idx < 4; ++idx) {
                out[idx] = observer_state.tEventList[idx];
                out[4 + idx] = observer_state.xEventList[idx];
            }
            out[8] = (realtype)observer_state.eventcount;
        }
        """
        ).replace("__EVENT_DIRECTION__", str(int(event_direction))),
        extra_options=("-DUSE_OBSERVER_THRESHOLD_CROSSING", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float32)
    x_values = np.array([-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0], dtype=np.float32)
    out = np.empty(9, dtype=np.float32)

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
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_threshold1_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    observed_times = out[: len(expected_times)]
    np.testing.assert_allclose(observed_times, expected_times, atol=1e-6, rtol=0.0)
    assert np.all(out[len(expected_times) : 4] == 0.0)
    np.testing.assert_allclose(
        out[4 : 4 + len(expected_times)],
        np.full(len(expected_times), 0.5, dtype=np.float32),
        atol=1e-6,
        rtol=0.0,
    )
    assert np.all(out[4 + len(expected_times) : 8] == 0.0)
    assert out[8] == pytest.approx(float(len(expected_times)))


@pytest.mark.parametrize(
    ("event_direction", "expected_times"),
    [
        (EventDirection.rising, np.array([1.5, 5.5], dtype=np.float32)),
        (EventDirection.falling, np.array([2.5], dtype=np.float32)),
        (EventDirection.either, np.array([1.5, 2.5, 5.5], dtype=np.float32)),
    ],
)
def test_normalized_threshold_crossing_observer_kernel_tracks_normalized_crossings(
    event_direction: EventDirection,
    expected_times: np.ndarray,
) -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        (
            """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
            0,
            0,
            8,
            __EVENT_DIRECTION__,
            RCONST(0.5),
            ZERO,
            ZERO,
            RCONST(0.75),
            ZERO,
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_normalized_threshold_crossing_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {ZERO};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            observer_state.xGlobalMax = ONE;
            observer_state.xGlobalMin = -ONE;
            initializeEventDetector(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            for (int idx = 0; idx < 4; ++idx)
                out[idx] = observer_state.tEventList[idx];
            out[4] = (realtype)observer_state.eventcount;
        }
        """
        ).replace("__EVENT_DIRECTION__", str(int(event_direction))),
        extra_options=("-DUSE_OBSERVER_NORMALIZED_THRESHOLD_CROSSING", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0], dtype=np.float32)
    x_values = np.array([-1.0, 0.0, 1.0, 0.0, -1.0, 0.0, 1.0], dtype=np.float32)
    out = np.empty(5, dtype=np.float32)

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
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_normalized_threshold_crossing_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    observed_times = out[: len(expected_times)]
    np.testing.assert_allclose(observed_times, expected_times, atol=1e-6, rtol=0.0)
    assert np.all(out[len(expected_times) : 4] == 0.0)
    assert out[4] == pytest.approx(float(len(expected_times)))


def test_local_maximum_observer_kernel_tracks_maxima_only() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        f"""
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {{
            0,
            0,
            4,
            0,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO,
            ZERO
        }};

        __kernel void run_local_maximum_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *dx_values,
            __global realtype *out,
            const uint n_values
        ) {{
            realtype ti = times[0];
            realtype xi[1] = {{x_values[0]}};
            realtype dxi[1] = {{dx_values[0]}};
            realtype auxi[1] = {{ZERO}};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {{
                ti = times[idx];
                xi[0] = x_values[idx];
                dxi[0] = dx_values[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {{
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }}
            }}

            out[0] = observer_state.tMaxList[0];
            out[1] = observer_state.xMaxList[0];
            out[2] = observer_state.tMaxList[1];
            out[3] = observer_state.xMaxList[1];
            out[4] = (realtype)observer_state.eventcount;
        }}
        """,
        extra_options=("-DUSE_OBSERVER_LOCAL_MAX", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    dt = 0.2
    sample_times = np.arange(0.0, 4.0 * math.pi + dt, dt, dtype=np.float64)
    x_values = np.sin(sample_times).astype(np.float32)
    dx_values = np.cos(sample_times).astype(np.float32)
    out = np.empty(5, dtype=np.float32)

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

    program.run_local_maximum_observer(
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
        np.array([math.pi / 2.0, 1.0, 5.0 * math.pi / 2.0, 1.0, 2.0], dtype=np.float32),
        atol=2e-2,
        rtol=0.0,
    )


def test_local_max_observer_kernel_uses_three_sample_extremum_helpers() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
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
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
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


def test_neighborhood_return_observer_kernel_interpolates_exit_events() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
            0,
            0,
            4,
            0,
            ZERO,
            ZERO,
            RCONST(0.25),
            ZERO,
            RCONST(0.25),
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_normalized_neighborhood_return_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[1] = {x_values[0]};
            realtype dxi[1] = {ZERO};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                warmupObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            }

            ti = times[0];
            xi[0] = x_values[0];
            initializeEventDetector(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tEventList[0];
            out[1] = observer_state.xEventList[0];
            out[2] = observer_state.tEventList[1];
            out[3] = observer_state.xEventList[1];
            out[4] = (realtype)observer_state.eventcount;
        }
        """,
        extra_options=("-DUSE_OBSERVER_NORMALIZED_NEIGHBORHOOD_RETURN", "-DN_VAR=1", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.arange(8, dtype=np.float32)
    x_values = np.array([0.0, -0.2, -0.6, -1.0, 0.0, -0.8, -1.0, 0.0], dtype=np.float32)
    out = np.empty(5, dtype=np.float32)

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
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_normalized_neighborhood_return_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(
        out,
        np.array([3.25, -0.75, 6.25, -0.75, 2.0], dtype=np.float32),
        atol=1e-6,
        rtol=0.0,
    )


def test_neighborhood_return_observer_interpolates_using_full_state_geometry() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"
        #include \"observers.cl\"

        __constant struct ObserverRuntimeSettings TEST_PARAMS = {
            0,
            0,
            4,
            0,
            ZERO,
            ZERO,
            RCONST(0.25),
            ZERO,
            RCONST(0.25),
            ZERO,
            ZERO,
            ZERO
        };

        __kernel void run_normalized_neighborhood_return_observer(
            __global const realtype *times,
            __global const realtype *x_values,
            __global const realtype *y_values,
            __global realtype *out,
            const uint n_values
        ) {
            realtype ti = times[0];
            realtype xi[2] = {x_values[0], y_values[0]};
            realtype dxi[2] = {ZERO, ZERO};
            realtype auxi[1] = {ZERO};
            ObserverState observer_state;

            initializeObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                xi[1] = y_values[idx];
                warmupObserverState(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
            }

            ti = times[0];
            xi[0] = x_values[0];
            xi[1] = y_values[0];
            initializeEventDetector(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);

            for (uint idx = 1; idx < n_values; ++idx) {
                ti = times[idx];
                xi[0] = x_values[idx];
                xi[1] = y_values[idx];
                updateObserverState(
                    &ti,
                    xi,
                    dxi,
                    auxi,
                    times[idx] - times[idx - 1],
                    times[idx] - times[0],
                    &observer_state,
                    &TEST_PARAMS
                );
                if (eventFunction(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS)) {
                    computeEventFeatures(&ti, xi, dxi, auxi, &observer_state, &TEST_PARAMS);
                }
            }

            out[0] = observer_state.tEventList[0];
            out[1] = (realtype)observer_state.eventcount;
        }
        """,
        extra_options=("-DUSE_OBSERVER_NORMALIZED_NEIGHBORHOOD_RETURN", "-DN_VAR=2", "-DN_AUX=0", "-DN_STORE_EVENTS=4"),
    )

    sample_times = np.arange(5, dtype=np.float32)
    x_values = np.array([0.0, -0.6, -1.0, -1.0, -1.0], dtype=np.float32)
    y_values = np.array([0.0, 0.0, 0.0, 1.0, 0.0], dtype=np.float32)
    out = np.empty(2, dtype=np.float32)

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
    y_buffer = pyopencl.Buffer(
        runtime.context,
        pyopencl.mem_flags.READ_ONLY | pyopencl.mem_flags.COPY_HOST_PTR,
        hostbuf=y_values,
    )
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_normalized_neighborhood_return_observer(
        runtime.queue,
        (1,),
        None,
        time_buffer,
        x_buffer,
        y_buffer,
        out_buffer,
        np.uint32(len(sample_times)),
    )
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    np.testing.assert_allclose(
        out,
        np.array([2.25, 1.0], dtype=np.float32),
        atol=1e-6,
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


def test_compensated_time_value_reconstructs_large_elapsed_with_kahan_sign() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void probe_compensated_time_value(__global realtype *out) {
            realtype elapsed = RCONST(4096.0);
            realtype elapsed_correction = ZERO;

            compensatedTimeAdd(&elapsed, &elapsed_correction, RCONST(2.0e-4));

            out[0] = elapsed;
            out[1] = elapsed_correction;
            out[2] = compensatedTimeValue(elapsed, elapsed_correction);
            out[3] = compensatedTimeFromOrigin(ZERO, elapsed, elapsed_correction);
            out[4] = compensatedTimeFromOriginAfterStep(ZERO, RCONST(4096.0), ZERO, RCONST(2.0e-4));
        }
        """,
    )

    out = np.empty(5, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.probe_compensated_time_value(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(np.float32(4096.0))
    assert out[1] == pytest.approx(np.float32(-2.0e-4), abs=1e-7)
    assert out[2] == pytest.approx(np.float32(4096.0))
    assert out[3] == pytest.approx(np.float32(4096.0))
    assert out[4] == pytest.approx(np.float32(4096.0))


def test_compensated_time_helpers_reduce_large_origin_failure_and_zero_origin_drift() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        __kernel void compare_time_prototypes(__global realtype *out) {
            const realtype dt = RCONST(0.01);

            realtype large_direct = RCONST(1000000.0);
            realtype large_elapsed = ZERO;
            realtype large_elapsed_correction = ZERO;
            for (uint step = 0; step < 10000; ++step) {
                large_direct += dt;
                compensatedTimeAdd(&large_elapsed, &large_elapsed_correction, dt);
            }
            out[0] = large_direct;
            out[1] = compensatedTimeFromOrigin(
                RCONST(1000000.0),
                large_elapsed,
                large_elapsed_correction
            );
            out[2] = compensatedTimeValue(large_elapsed, large_elapsed_correction);
            out[3] = compensatedTimeFromOriginAfterStep(
                RCONST(1000000.0),
                large_elapsed,
                large_elapsed_correction,
                ZERO
            );

            realtype zero_direct = ZERO;
            realtype zero_elapsed = ZERO;
            realtype zero_elapsed_correction = ZERO;
            for (uint step = 0; step < 200000; ++step) {
                zero_direct += dt;
                compensatedTimeAdd(&zero_elapsed, &zero_elapsed_correction, dt);
            }
            out[4] = zero_direct;
            out[5] = compensatedTimeValue(zero_elapsed, zero_elapsed_correction);
            out[6] = compensatedTimeFromOrigin(ZERO, zero_elapsed, zero_elapsed_correction);
            out[7] = compensatedTimeFromOriginAfterStep(
                ZERO,
                zero_elapsed,
                zero_elapsed_correction,
                ZERO
            );

            realtype step3_elapsed = ZERO;
            realtype step3_elapsed_correction = ZERO;
            for (uint step = 0; step < 3; ++step)
                compensatedTimeAdd(&step3_elapsed, &step3_elapsed_correction, dt);

            out[8] = compensatedTimeFromOrigin(
                RCONST(1000000.0),
                step3_elapsed,
                step3_elapsed_correction
            );
            out[9] = compensatedTimeFromOriginAfterStep(
                RCONST(1000000.0),
                step3_elapsed,
                step3_elapsed_correction,
                RCONST(0.5) * dt
            );
        }
        """,
    )

    out = np.empty(10, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.compare_time_prototypes(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(np.float32(1.0e6))
    assert out[1] == pytest.approx(np.float32(1000100.0))
    assert out[2] == pytest.approx(np.float32(100.0))
    assert out[3] == pytest.approx(np.float32(1000100.0))
    assert abs(float(out[4]) - 2000.0) > 1.0
    assert out[5] == pytest.approx(np.float32(2000.0))
    assert out[6] == pytest.approx(np.float32(2000.0))
    assert out[7] == pytest.approx(np.float32(2000.0))
    assert out[8] == pytest.approx(np.float32(1000000.0))
    assert out[9] == pytest.approx(np.float32(1000000.0625))
    assert out[9] > out[8]


def test_adaptive_he12_component_step_uses_compensated_relative_time_at_large_origin() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_utilities.cl\"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        );

        #include \"steppers/adaptive_he12.clh\"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        ) {
            (void)t;
            (void)x_;
            (void)p_;
            (void)aux_;
            (void)w_;
            dx_[0] = ONE;
        }

        __kernel void compare_large_origin_adaptive_step(__global realtype *out) {
            const realtype t_origin = RCONST(1000000.0);
            const realtype dt = RCONST(0.01);

            realtype naive_t = t_origin;
            realtype compensated_t = t_origin;
            realtype solve_elapsed = ZERO;
            realtype solve_elapsed_correction = ZERO;

            realtype naive_x[1] = {ZERO};
            realtype compensated_x[1] = {ZERO};
            realtype naive_k1[1] = {ONE};
            realtype compensated_k1[1] = {ONE};
            realtype pars[1] = {ZERO};
            realtype aux[1] = {ZERO};
            realtype err[1] = {ZERO};
            realtype wi[1] = {ZERO};

            for (uint step = 0; step < 10000; ++step) {
                realtype t_new = naive_t + dt;
                realtype effective_dt = t_new - naive_t;
                naive_x[0] += effective_dt;
                naive_t = t_new;

                do_step(
                    t_origin,
                    &solve_elapsed,
                    &solve_elapsed_correction,
                    &compensated_t,
                    compensated_x,
                    compensated_k1,
                    pars,
                    dt,
                    aux,
                    err,
                    wi
                );
            }

            out[0] = naive_t;
            out[1] = compensated_t;
            out[2] = compensatedTimeValue(solve_elapsed, solve_elapsed_correction);
            out[3] = naive_x[0];
            out[4] = compensated_x[0];
        }
        """,
        extra_options=("-DN_VAR=1", "-DN_AUX=0", "-DN_PAR=1", "-DN_WIENER=0"),
    )

    out = np.empty(5, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.compare_large_origin_adaptive_step(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(np.float32(1.0e6))
    assert out[1] == pytest.approx(np.float32(1000100.0))
    assert out[2] == pytest.approx(np.float32(100.0), abs=1e-4)
    assert out[3] == pytest.approx(np.float32(0.0))
    assert out[4] == pytest.approx(np.float32(100.0), abs=5e-3)


def test_fixed_stepper_reports_accepted_step_width_without_elapsed_differencing() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include \"clODE_random.cl\"
        #include \"clODE_struct_defs.cl\"
        #include \"clODE_utilities.cl\"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        );

        #include \"steppers.cl\"

        __constant struct IntegrationSettings TEST_SETTINGS = {
            RCONST(1.0e-4),
            RCONST(1.0e-4),
            ZERO,
            ZERO,
            8
        };
        __constant realtype TEST_TSPAN[2] = {RCONST(1000000.0), RCONST(1005000.0)};

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        ) {
            (void)t;
            (void)x_;
            (void)p_;
            (void)aux_;
            (void)w_;
            dx_[0] = ONE;
        }

        __kernel void probe_fixed_stepper_width(__global realtype *out) {
            realtype ti = RCONST(1004096.0);
            realtype dt = RCONST(1.0e-4);
            realtype accepted_dt = ZERO;
            realtype solve_elapsed = RCONST(4096.0);
            realtype solve_elapsed_correction = ZERO;
            realtype xi[1] = {ZERO};
            realtype k1[1] = {ONE};
            realtype pars[1] = {ZERO};
            realtype aux[1] = {ZERO};
            realtype wi[1] = {ZERO};
            struct rngData rd;

            for (int j = 0; j < N_RNGSTATE; ++j)
                rd.state[j] = 0UL;
            rd.randnUselast = 0;
            rd.randnLast = ZERO;

            realtype previous_elapsed = compensatedTimeValue(
                solve_elapsed,
                solve_elapsed_correction
            );
            int stepflag = stepper(
                &ti,
                &solve_elapsed,
                &solve_elapsed_correction,
                xi,
                k1,
                pars,
                &TEST_SETTINGS,
                &dt,
                &accepted_dt,
                TEST_TSPAN[0],
                TEST_TSPAN[1] - TEST_TSPAN[0],
                aux,
                wi,
                &rd
            );

            out[0] = stepflag;
            out[1] = accepted_dt;
            out[2] = compensatedTimeValue(solve_elapsed, solve_elapsed_correction) - previous_elapsed;
            out[3] = xi[0];
            out[4] = ti;
        }
        """,
        extra_options=("-DEXPLICIT_EULER", "-DN_VAR=1", "-DN_AUX=0", "-DN_PAR=1", "-DN_WIENER=0"),
    )

    out = np.empty(5, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.probe_fixed_stepper_width(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    assert out[0] == pytest.approx(np.float32(4.0))
    assert out[1] == pytest.approx(np.float32(0.0))
    assert out[2] == pytest.approx(np.float32(0.0))
    assert out[3] == pytest.approx(np.float32(0.0))
    assert out[4] == pytest.approx(np.float32(1004096.0))


def test_fixed_heun_component_step_uses_origin_plus_elapsed_end_stage_time() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include "clODE_random.cl"
        #include "clODE_utilities.cl"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        );

        #include "steppers/fixed_explicit_Trapezoidal.clh"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        ) {
            (void)x_;
            (void)aux_;
            (void)w_;
            dx_[0] = p_[0] * (t - p_[1]);
        }

        __kernel void run_fixed_heun_step(__global realtype *out) {
            const realtype t_origin = ONE;
            realtype elapsed = ONE;
            realtype elapsed_correction = ZERO;
            const realtype dt = RCONST(1.0e-4);
            realtype ti = compensatedTimeFromOrigin(t_origin, elapsed, elapsed_correction);
            realtype xi[1] = {ZERO};
            realtype pars[2] = {RCONST(1.0e8), RCONST(2.0)};
            realtype aux[1] = {ZERO};
            realtype wi[1] = {ZERO};
            realtype k1[1] = {ZERO};
            struct rngData rd;

            for (int j = 0; j < N_RNGSTATE; ++j)
                rd.state[j] = 0UL;
            rd.randnUselast = 0;
            rd.randnLast = ZERO;

            getRHS(
                ti,
                xi,
                pars,
                k1,
                aux,
                wi
            );
            do_step(
                t_origin,
                &elapsed,
                &elapsed_correction,
                &ti,
                xi,
                k1,
                pars,
                dt,
                aux,
                wi,
                &rd
            );
            out[0] = xi[0];
            out[1] = k1[0];
            out[2] = ti;
        }
        """,
        extra_options=("-DN_VAR=1", "-DN_AUX=0", "-DN_PAR=2", "-DN_WIENER=0"),
    )

    out = np.empty(3, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_fixed_heun_step(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    origin = np.float32(1.0)
    elapsed = np.float32(1.0)
    correction = np.float32(0.0)
    dt = np.float32(1.0e-4)
    scale = np.float32(1.0e8)
    anchor = np.float32(2.0)
    t_start = _reference_compensated_time_from_origin(origin, elapsed, correction)
    t_end = _reference_compensated_time_from_origin_after_step(origin, elapsed, correction, dt)
    expected = np.float32(
        dt * np.float32(0.5) * (scale * (t_start - anchor) + scale * (t_end - anchor))
    )

    assert out[0] == pytest.approx(expected, abs=1e-6)
    assert out[1] == pytest.approx(scale * (t_end - anchor), abs=1e-3)
    assert out[2] == pytest.approx(t_end, abs=1e-7)


def test_fixed_rk4_component_step_uses_origin_plus_elapsed_stage_times() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    program = _build_program(
        runtime,
        """
        #include "clODE_random.cl"
        #include "clODE_utilities.cl"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        );

        #include "steppers/fixed_explicit_RK4.clh"

        void getRHS(
            const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]
        ) {
            (void)x_;
            (void)aux_;
            (void)w_;
            dx_[0] = p_[0] * (t - p_[1]);
        }

        __kernel void run_fixed_rk4_step(__global realtype *out) {
            const realtype t_origin = RCONST(100.0);
            realtype elapsed = RCONST(10.0);
            realtype elapsed_correction = ZERO;
            const realtype dt = RCONST(1.0e-4);
            realtype ti = compensatedTimeFromOrigin(t_origin, elapsed, elapsed_correction);
            realtype xi[1] = {ZERO};
            realtype pars[2] = {RCONST(1.0e8), RCONST(110.0)};
            realtype aux[1] = {ZERO};
            realtype wi[1] = {ZERO};
            realtype k1[1] = {ZERO};
            struct rngData rd;

            for (int j = 0; j < N_RNGSTATE; ++j)
                rd.state[j] = 0UL;
            rd.randnUselast = 0;
            rd.randnLast = ZERO;

            getRHS(
                ti,
                xi,
                pars,
                k1,
                aux,
                wi
            );
            do_step(
                t_origin,
                &elapsed,
                &elapsed_correction,
                &ti,
                xi,
                k1,
                pars,
                dt,
                aux,
                wi,
                &rd
            );
            out[0] = xi[0];
            out[1] = k1[0];
            out[2] = ti;
        }
        """,
        extra_options=("-DN_VAR=1", "-DN_AUX=0", "-DN_PAR=2", "-DN_WIENER=0"),
    )

    out = np.empty(3, dtype=np.float32)
    out_buffer = pyopencl.Buffer(runtime.context, pyopencl.mem_flags.WRITE_ONLY, out.nbytes)

    program.run_fixed_rk4_step(runtime.queue, (1,), None, out_buffer)
    pyopencl.enqueue_copy(runtime.queue, out, out_buffer).wait()

    origin = np.float32(100.0)
    elapsed = np.float32(10.0)
    correction = np.float32(0.0)
    dt = np.float32(1.0e-4)
    scale = np.float32(1.0e8)
    anchor = np.float32(110.0)
    t_start = _reference_compensated_time_from_origin(origin, elapsed, correction)
    t_mid = _reference_compensated_time_from_origin_after_step(
        origin,
        elapsed,
        correction,
        np.float32(0.5) * dt,
    )
    t_end = _reference_compensated_time_from_origin_after_step(origin, elapsed, correction, dt)
    expected = np.float32(
        dt
        * (
            scale * (t_start - anchor)
            + np.float32(2.0) * scale * (t_mid - anchor)
            + np.float32(2.0) * scale * (t_mid - anchor)
            + scale * (t_end - anchor)
        )
        / np.float32(6.0)
    )

    assert out[0] == pytest.approx(expected, abs=1e-6)
    assert out[1] == pytest.approx(scale * (t_end - anchor), abs=1e-3)
    assert out[2] == pytest.approx(t_end, abs=1e-7)