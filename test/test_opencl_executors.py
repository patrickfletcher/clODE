import numpy as np
import pytest

pytest.importorskip("pyopencl")

import clode

from clode._opencl import (
    OpenCLFeatureExecutor,
    OpenCLRuntime,
    OpenCLTrajectoryExecutor,
    OpenCLTransientExecutor,
)
from clode.observers.types import EventOutputSettings, ObserverRuntimeSettings
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


def _make_feature_executor() -> OpenCLFeatureExecutor:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    executor = OpenCLFeatureExecutor(
        ProblemInfo("stable_linear.cl", ["x", "y"], ["a", "b"], [], 0),
        load_rhs_source(model_path("stable_linear.cl")),
        "rk4",
        "summary",
        ObserverRuntimeSettings(f_var_ix=0),
        EventOutputSettings(),
        None,
        True,
        runtime,
        _clode_root_dir,
    )
    executor.build_cl()
    executor.set_solver_params(
        SolverParams(
            dt=REQUESTED_DT,
            dtmax=REQUESTED_DT,
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


def _make_fixed_transient_executor(*, max_steps: int) -> OpenCLTransientExecutor:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    executor = OpenCLTransientExecutor(
        ProblemInfo("stable_linear.cl", ["x", "y"], ["a", "b"], [], 0),
        load_rhs_source(model_path("stable_linear.cl")),
        "rk4",
        True,
        runtime,
        _clode_root_dir,
    )
    executor.build_cl()
    executor.set_solver_params(
        SolverParams(
            dt=0.3,
            dtmax=0.3,
            abstol=1e-7,
            reltol=1e-6,
            max_steps=max_steps,
            max_store=1,
            nout=1,
        )
    )
    executor.set_tspan((0.0, 1.0))
    executor.set_problem_data([2.0, -1.5], [0.5, 1.25])
    return executor


def test_problem_data_reset_restores_requested_dt() -> None:
    executor = _make_transient_executor()

    executor.transient()
    continued_dt = np.asarray(executor.get_dt(), dtype=np.float64)
    accepted_dt = np.asarray(executor.get_last_accepted_dt(), dtype=np.float64)

    assert continued_dt[0] < REQUESTED_DT / 10.0
    assert accepted_dt[0] > continued_dt[0]
    assert executor.get_status() == [clode.SolverStatus.COMPLETED]

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


def test_trajectory_stride_change_preserves_solver_state_and_reuses_trajectory_buffers() -> None:
    executor = _make_trajectory_executor()

    executor.trajectory()
    continued_dt = np.asarray(executor.get_dt(), dtype=np.float64)
    step_count = np.asarray(executor.get_step_count(), dtype=np.uint64)
    accepted_dt = np.asarray(executor.get_last_accepted_dt(), dtype=np.float64)
    final_time = np.asarray(executor.get_tf(), dtype=np.float64)
    trajectory_buffers = executor._trajectory_buffers

    assert continued_dt[0] > 0.02
    assert step_count[0] > 0
    assert accepted_dt[0] > 0.0
    assert executor.get_status() == [clode.SolverStatus.OUTPUT_CAPACITY_REACHED]
    assert trajectory_buffers is not None

    executor.set_solver_params(
        SolverParams(
            dt=0.01,
            dtmax=0.2,
            abstol=1e-9,
            reltol=1e-8,
            max_steps=1024,
            max_store=16,
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
    assert executor._trajectory_buffers is trajectory_buffers
    assert executor.get_t() == []


def test_trajectory_capacity_change_preserves_solver_state_and_invalidates_trajectory_buffers() -> None:
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
    assert executor._trajectory_buffers is None
    assert executor.get_t() == []


def test_feature_runtime_setting_change_preserves_program_and_buffers() -> None:
    executor = _make_feature_executor()

    executor.features()
    feature_buffers = executor._feature_buffers
    program_bundle = executor._program_bundle

    assert feature_buffers is not None
    assert program_bundle is not None
    assert executor.get_status() == [clode.SolverStatus.COMPLETED]
    assert executor.get_step_count()[0] > 0
    assert executor.get_last_accepted_dt()[0] > 0.0
    assert executor.get_feature_names()[0] == "max x"

    updated_params = executor.get_observer_params()
    updated_params.f_var_ix = 1
    executor.set_observer_params(updated_params)

    assert executor._feature_buffers is feature_buffers
    assert executor._program_bundle is program_bundle
    assert executor.get_feature_names()[0] == "max x"
    assert executor.get_f() == []


def test_fixed_step_transient_status_reports_max_steps_exhaustion() -> None:
    executor = _make_fixed_transient_executor(max_steps=1)

    executor.transient()

    assert executor.get_status() == [clode.SolverStatus.MAX_STEPS_REACHED]
    assert executor.get_step_count() == [1]
    np.testing.assert_allclose(
        np.asarray(executor.get_last_accepted_dt(), dtype=np.float64),
        [0.3],
        atol=1e-7,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        np.asarray(executor.get_tf(), dtype=np.float64),
        [0.3],
        atol=1e-7,
        rtol=0.0,
    )