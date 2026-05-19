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