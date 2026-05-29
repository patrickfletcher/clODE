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
FIXED_ATOL = 1e-6
ADAPTIVE_ATOL = 1e-4
CONTINUATION_DT = 0.125
CONTINUATION_ATOL = 1e-12


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


def _feature_matrix(output: clode.ObserverOutput) -> np.ndarray:
    return np.asarray(output.to_ndarray(), dtype=np.float64)


def test_rk4_summary_linear_state_statistics_match_exact_values() -> None:
    simulator = _settled_feature_simulator("stable_linear")
    output = simulator.features()

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
    assert int(simulator.get_step_count().reshape(-1)[0]) == fixed_step_step_count(
        FIXED_SETTLE_END,
        FIXED_SETTLE_END + FIXED_WINDOW,
        FIXED_DT,
    )


def test_rk4_summary_linear_aux_statistics_match_exact_values() -> None:
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


def test_dormand_prince_summary_hopf_cycle_statistics_match_exact_values() -> None:
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


def test_dormand_prince_summary_large_origin_matches_zero_origin_statistics() -> None:
    duration = 4.0 * hopf_cycle_period()
    kwargs = dict(
        model_name="hopf_normal_form",
        stepper=clode.Stepper.dormand_prince,
        dt=0.05,
        dtmax=0.1,
        abstol=1e-7,
        reltol=1e-6,
        max_steps=4096,
        single_precision=True,
    )

    zero_origin = make_feature_simulator(
        t_span=(0.0, duration),
        **kwargs,
    ).features()
    large_origin_simulator = make_feature_simulator(
        t_span=(10000.0, 10000.0 + duration),
        **kwargs,
    )
    large_origin = large_origin_simulator.features()

    assert zero_origin is not None
    assert large_origin is not None

    np.testing.assert_allclose(
        _feature_matrix(large_origin),
        _feature_matrix(zero_origin),
        atol=3e-4,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        large_origin_simulator.get_final_time(),
        [10000.0 + duration],
        atol=3e-4,
        rtol=0.0,
    )


def test_rk4_summary_continuation_matches_single_run() -> None:
    full = make_feature_simulator(
        "stable_linear_aux",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=CONTINUATION_DT,
        dtmax=CONTINUATION_DT,
        max_steps=64,
        single_precision=False,
    )
    split = make_feature_simulator(
        "stable_linear_aux",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 0.5),
        dt=CONTINUATION_DT,
        dtmax=CONTINUATION_DT,
        max_steps=64,
        single_precision=False,
    )

    full_output = full.features()
    split.features()
    first_final_time = float(split.get_final_time()[0])
    split.set_tspan((first_final_time, 1.0))
    split_output = split.features()

    assert full_output is not None
    assert split_output is not None

    np.testing.assert_allclose(
        _feature_matrix(split_output),
        _feature_matrix(full_output),
        atol=CONTINUATION_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        split.get_final_time(),
        full.get_final_time(),
        atol=CONTINUATION_ATOL,
        rtol=0.0,
    )
