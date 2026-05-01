from __future__ import annotations

import math

import numpy as np

import clode
from test.core_numerics.helpers import make_trajectory_simulator, structured_to_array
from test.core_numerics.reference import (
    fixed_step_stored_times,
    hopf_on_cycle_aux,
    hopf_on_cycle_derivative,
    hopf_on_cycle_state,
    stable_linear_aux,
    stable_linear_derivative,
    stable_linear_state,
)

FIXED_DT = 0.05
FIXED_MAX_STEPS = 64
FIXED_ATOL = 1e-5
ADAPTIVE_ATOL = 1e-5


def test_rk4_stable_linear_trajectory_matches_exact_samples() -> None:
    trajectory = make_trajectory_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=FIXED_MAX_STEPS,
        nout=1,
    ).trajectory()

    times = np.asarray(trajectory.t, dtype=np.float64)
    expected_state = stable_linear_state(times)
    expected_derivative = stable_linear_derivative(times)

    np.testing.assert_allclose(
        structured_to_array(trajectory, "x"),
        expected_state,
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        structured_to_array(trajectory, "dx"),
        expected_derivative,
        atol=FIXED_ATOL,
        rtol=0.0,
    )


def test_rk4_linear_aux_trajectory_returns_exact_aux_samples() -> None:
    trajectory = make_trajectory_simulator(
        "stable_linear_aux",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=FIXED_MAX_STEPS,
        nout=1,
    ).trajectory()

    times = np.asarray(trajectory.t, dtype=np.float64)
    expected_aux = stable_linear_aux(times)

    np.testing.assert_allclose(structured_to_array(trajectory, "aux"), expected_aux, atol=FIXED_ATOL, rtol=0.0)


def test_trajectory_nout_and_max_store_contract() -> None:
    simulator = make_trajectory_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.25,
        dtmax=0.25,
        max_steps=32,
        max_store=3,
        nout=2,
    )
    trajectory = simulator.trajectory()

    expected_times = fixed_step_stored_times(0.0, 1.0, 0.25, 2)

    np.testing.assert_allclose(np.asarray(trajectory.t, dtype=np.float64), expected_times, atol=0.0, rtol=0.0)
    assert int(simulator._device_n_stored[0]) == len(expected_times) - 1


def test_dormand_prince_hopf_trajectory_matches_exact_cycle() -> None:
    trajectory = make_trajectory_simulator(
        "hopf_normal_form",
        stepper=clode.Stepper.dormand_prince,
        t_span=(0.0, 2.0 * math.pi),
        dt=0.05,
        dtmax=0.1,
        abstol=1e-7,
        reltol=1e-6,
        max_steps=2048,
        max_store=2048,
        nout=1,
    ).trajectory()

    times = np.asarray(trajectory.t, dtype=np.float64)
    expected_state = hopf_on_cycle_state(times)
    expected_derivative = hopf_on_cycle_derivative(times)
    expected_aux = hopf_on_cycle_aux(times)

    np.testing.assert_allclose(
        structured_to_array(trajectory, "x"),
        expected_state,
        atol=ADAPTIVE_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        structured_to_array(trajectory, "dx"),
        expected_derivative,
        atol=ADAPTIVE_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        structured_to_array(trajectory, "aux"),
        expected_aux,
        atol=ADAPTIVE_ATOL,
        rtol=0.0,
    )
