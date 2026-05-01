from __future__ import annotations

import numpy as np

import clode
from test.core_numerics.helpers import make_feature_simulator
from test.core_numerics.reference import (
    fixed_step_observer_times,
    fixed_step_step_count,
    hopf_cycle_period,
    stable_linear_aux,
    stable_linear_derivative,
    stable_linear_state,
)

FIXED_DT = 0.05
FIXED_SETTLE_END = 8.0
FIXED_WINDOW = 1.0
FIXED_MAX_STEPS = 256
FIXED_ATOL = 1e-5
ADAPTIVE_ATOL = 1e-4


def _settled_feature_simulator(model_name: str) -> clode.FeatureSimulator:
    simulator = make_feature_simulator(
        model_name,
        stepper=clode.Stepper.rk4,
        t_span=(0.0, FIXED_SETTLE_END),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )
    simulator.transient(update_x0=True)
    simulator.set_tspan((FIXED_SETTLE_END, FIXED_SETTLE_END + FIXED_WINDOW))
    return simulator


def test_rk4_basicall_linear_state_statistics_match_exact_values() -> None:
    output = _settled_feature_simulator("stable_linear").features()

    assert output is not None
    sample_times = fixed_step_observer_times(
        FIXED_SETTLE_END,
        FIXED_SETTLE_END + FIXED_WINDOW,
        FIXED_DT,
    )
    states = stable_linear_state(sample_times)
    derivatives = stable_linear_derivative(sample_times)

    np.testing.assert_allclose(
        output.get_var_max("x"),
        states[:, 0].max(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_min("x"),
        states[:, 0].min(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_mean("x"),
        states[:, 0].mean(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_max("y"),
        states[:, 1].max(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_min("y"),
        states[:, 1].min(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_mean("y"),
        states[:, 1].mean(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_max_slope("x"),
        derivatives[:, 0].max(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_min_slope("x"),
        derivatives[:, 0].min(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_max_slope("y"),
        derivatives[:, 1].max(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_min_slope("y"),
        derivatives[:, 1].min(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    assert output.get_var_count("step") == fixed_step_step_count(
        FIXED_SETTLE_END,
        FIXED_SETTLE_END + FIXED_WINDOW,
        FIXED_DT,
    )


def test_rk4_basicall_linear_aux_statistics_match_exact_values() -> None:
    output = _settled_feature_simulator("stable_linear_aux").features()

    assert output is not None
    sample_times = fixed_step_observer_times(
        FIXED_SETTLE_END,
        FIXED_SETTLE_END + FIXED_WINDOW,
        FIXED_DT,
    )
    aux = stable_linear_aux(sample_times)

    np.testing.assert_allclose(
        output.get_var_max("sum"),
        aux[:, 0].max(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_min("sum"),
        aux[:, 0].min(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_mean("sum"),
        aux[:, 0].mean(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_max("combo"),
        aux[:, 1].max(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_min("combo"),
        aux[:, 1].min(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        output.get_var_mean("combo"),
        aux[:, 1].mean(),
        atol=FIXED_ATOL,
        rtol=0.0,
    )


def test_dormand_prince_basicall_hopf_cycle_statistics_match_exact_values() -> None:
    output = make_feature_simulator(
        "hopf_normal_form",
        stepper=clode.Stepper.dormand_prince,
        t_span=(0.0, 4.0 * hopf_cycle_period()),
        dt=0.05,
        dtmax=0.1,
        abstol=1e-7,
        reltol=1e-6,
        max_steps=4096,
    ).features()

    assert output is not None

    np.testing.assert_allclose(output.get_var_max("x"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_min("x"), -1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_mean("x"), 0.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_max("y"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_min("y"), -1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_mean("y"), 0.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_max_slope("x"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_min_slope("x"), -1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_max_slope("y"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_min_slope("y"), -1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_max("r2"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_min("r2"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
    np.testing.assert_allclose(output.get_var_mean("r2"), 1.0, atol=ADAPTIVE_ATOL, rtol=0.0)
