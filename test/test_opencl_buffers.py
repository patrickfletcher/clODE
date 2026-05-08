import numpy as np
import pytest

pytest.importorskip("pyopencl")

from clode._opencl import ArrayLayout, BufferManager, OpenCLRuntime, Precision, ProblemShape
from clode.simulation import SolverParams
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def test_array_layout_matches_current_fortran_order_contract() -> None:
    matrix = np.array([[1.0, 2.0], [3.0, 4.0]])

    flattened = ArrayLayout.flatten_problem_matrix(matrix)
    reshaped = ArrayLayout.reshape_state(flattened, ensemble_size=2, width=2)

    assert flattened.tolist() == [1.0, 3.0, 2.0, 4.0]
    assert np.array_equal(reshaped, matrix)


def test_buffer_manager_allocates_common_buffers_with_current_sizes() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    manager = BufferManager(runtime, Precision.SINGLE)
    shape = ProblemShape(n_var=2, n_par=1, n_aux=0, n_wiener=1)

    buffers = manager.allocate_common(ensemble_size=3, shape=shape)

    assert buffers.tspan.size == 2 * np.dtype(np.float32).itemsize
    assert buffers.x0.size == 3 * 2 * np.dtype(np.float32).itemsize
    assert buffers.pars.size == 3 * 1 * np.dtype(np.float32).itemsize
    assert buffers.xf.size == 3 * 2 * np.dtype(np.float32).itemsize
    assert buffers.rng_state.size == 3 * 2 * np.dtype(np.uint64).itemsize
    assert buffers.dt.size == 3 * np.dtype(np.float32).itemsize
    assert buffers.tf.size == 3 * np.dtype(np.float32).itemsize


def test_buffer_manager_roundtrips_problem_data_and_dt() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    manager = BufferManager(runtime, Precision.DOUBLE)
    shape = ProblemShape(n_var=2, n_par=2, n_aux=0, n_wiener=0)
    buffers = manager.allocate_common(ensemble_size=3, shape=shape)

    initial_state = np.array(
        [[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]], dtype=np.float64
    )
    parameters = np.array(
        [[0.1, 1.1], [0.2, 1.2], [0.3, 1.3]], dtype=np.float64
    )
    dt_values = np.array([0.01, 0.02, 0.03], dtype=np.float64)
    rng_state = np.array(
        [[11, 21], [12, 22], [13, 23]],
        dtype=np.uint64,
    )

    manager.upload_problem_data(buffers, initial_state, parameters)
    manager.upload_dt(buffers, dt_values)
    manager.upload_rng_state(buffers, rng_state)

    assert np.array_equal(manager.download_x0(buffers), initial_state)
    assert np.array_equal(manager.download_pars(buffers), parameters)
    assert np.array_equal(manager.download_dt(buffers), dt_values)
    assert np.array_equal(manager.download_rng_state(buffers), rng_state)


def test_buffer_manager_packs_solver_params_with_current_struct_layout() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    manager = BufferManager(runtime, Precision.SINGLE)
    buffers = manager.allocate_common(
        ensemble_size=1,
        shape=ProblemShape(n_var=1, n_par=0, n_aux=0, n_wiener=0),
    )
    solver_params = SolverParams(
        dt=0.125,
        dtmax=0.5,
        abstol=1e-6,
        reltol=1e-3,
        max_steps=123,
        max_store=456,
        nout=7,
    )

    packed = manager.upload_solver_params(buffers, solver_params)

    assert packed.dtype.fields is not None
    assert packed.dtype.itemsize == 28
    assert float(packed["dt"]) == pytest.approx(0.125)
    assert int(packed["max_steps"]) == 123