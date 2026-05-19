import numpy as np
import pytest

pytest.importorskip("pyopencl")

from clode._opencl import (
    OpenCLRuntime,
    OpenCLTrajectoryExecutor,
    OpenCLTransientExecutor,
)
from clode.problem._core import ProblemInfo, load_rhs_source
from clode.runtime import _clode_root_dir
from clode.simulation import SolverParams
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID, model_path


REQUESTED_DT = 0.05


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def _make_transient_executor() -> OpenCLTransientExecutor:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    executor = OpenCLTransientExecutor(
        ProblemInfo("stable_linear.cl", ["x", "y"], ["a", "b"], [], 0),
        load_rhs_source(model_path("stable_linear.cl")),
        "dopri5",
        True,
        runtime,
        _clode_root_dir,
    )
    executor.build_cl()
    executor.set_solver_params(
        SolverParams(
            dt=REQUESTED_DT,
            dtmax=0.1,
            abstol=1e-7,
            reltol=1e-6,
            max_steps=256,
            max_store=1,
            nout=1,
        )
    )
    executor.set_tspan((0.0, 1.0))
    executor.set_problem_data([2.0, -1.5], [0.5, 1.25])
    return executor


def _make_trajectory_executor() -> OpenCLTrajectoryExecutor:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    executor = OpenCLTrajectoryExecutor(
        ProblemInfo("stable_linear.cl", ["x", "y"], ["a", "b"], [], 0),
        load_rhs_source(model_path("stable_linear.cl")),
        "dopri5",
        True,
        runtime,
        _clode_root_dir,
    )
    executor.build_cl()
    executor.set_solver_params(
        SolverParams(
            dt=0.01,
            dtmax=0.2,
            abstol=1e-9,
            reltol=1e-8,
            max_steps=1024,
            max_store=16,
            nout=1,
        )
    )
    executor.set_tspan((0.0, 4.0))
    executor.set_problem_data([2.0, -1.5], [0.5, 1.25])
    return executor


def test_problem_data_reset_restores_requested_dt() -> None:
    executor = _make_transient_executor()

    executor.transient()
    continued_dt = np.asarray(executor.get_dt(), dtype=np.float64)

    assert continued_dt[0] < REQUESTED_DT / 10.0

    executor.set_problem_data([2.0, -1.5], [0.5, 1.25])

    np.testing.assert_allclose(
        np.asarray(executor.get_dt(), dtype=np.float64),
        [REQUESTED_DT],
        atol=1e-7,
        rtol=0.0,
    )


def test_tspan_change_preserves_current_dt_and_invalidates_old_results() -> None:
    executor = _make_transient_executor()

    executor.transient()
    continued_dt = np.asarray(executor.get_dt(), dtype=np.float64)

    executor.set_tspan((1.0, 2.0))

    np.testing.assert_allclose(
        np.asarray(executor.get_dt(), dtype=np.float64),
        continued_dt,
        atol=0.0,
        rtol=0.0,
    )
    assert executor.get_xf() is None
    assert executor.get_tf() == []


def test_trajectory_output_change_preserves_solver_state_and_invalidates_trajectory_buffers() -> None:
    executor = _make_trajectory_executor()

    executor.trajectory()
    continued_dt = np.asarray(executor.get_dt(), dtype=np.float64)
    final_time = np.asarray(executor.get_tf(), dtype=np.float64)

    assert continued_dt[0] > 0.02

    executor.set_solver_params(
        SolverParams(
            dt=0.01,
            dtmax=0.2,
            abstol=1e-9,
            reltol=1e-8,
            max_steps=1024,
            max_store=4,
            nout=2,
        )
    )

    np.testing.assert_allclose(
        np.asarray(executor.get_dt(), dtype=np.float64),
        continued_dt,
        atol=0.0,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(executor.get_tf(), dtype=np.float64),
        final_time,
        atol=0.0,
        rtol=0.0,
    )

    with pytest.raises(RuntimeError, match="Trajectory buffers"):
        executor.get_t()