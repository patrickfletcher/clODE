from __future__ import annotations

import numpy as np

import clode
from test.core_numerics.helpers import make_simulator, set_repeat_ensemble_and_seed
from test.core_numerics.reference import ou_stationary_mean, ou_stationary_variance

REPRO_ENSEMBLE = 128
MOMENT_ENSEMBLE = 2048
REPRO_DT = 0.125
CONTINUATION_ATOL = 1e-12


def test_stochastic_euler_ou_seeded_ensemble_is_reproducible() -> None:
    first = make_simulator(
        "ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=REPRO_DT,
        dtmax=REPRO_DT,
        max_steps=64,
    )
    second = make_simulator(
        "ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=REPRO_DT,
        dtmax=REPRO_DT,
        max_steps=64,
    )

    set_repeat_ensemble_and_seed(first, ensemble_size=REPRO_ENSEMBLE, seed=123)
    set_repeat_ensemble_and_seed(second, ensemble_size=REPRO_ENSEMBLE, seed=123)

    first.transient()
    second.transient()

    np.testing.assert_array_equal(first.get_final_state(), second.get_final_state())
    np.testing.assert_array_equal(first.get_dt(), second.get_dt())
    np.testing.assert_array_equal(first.get_final_time(), second.get_final_time())


def test_stochastic_euler_ou_continuation_matches_seeded_single_run() -> None:
    full = make_simulator(
        "ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=REPRO_DT,
        dtmax=REPRO_DT,
        max_steps=64,
        single_precision=False,
    )
    split = make_simulator(
        "ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 0.5),
        dt=REPRO_DT,
        dtmax=REPRO_DT,
        max_steps=64,
        single_precision=False,
    )

    set_repeat_ensemble_and_seed(full, ensemble_size=REPRO_ENSEMBLE, seed=321)
    set_repeat_ensemble_and_seed(split, ensemble_size=REPRO_ENSEMBLE, seed=321)

    full.transient()
    split.transient()
    first_final_time = float(split.get_final_time().reshape(-1)[0])
    split.set_tspan((first_final_time, 1.0))
    split.transient()

    np.testing.assert_allclose(
        split.get_final_state(),
        full.get_final_state(),
        atol=CONTINUATION_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        split.get_dt(),
        full.get_dt(),
        atol=CONTINUATION_ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        split.get_final_time(),
        full.get_final_time(),
        atol=CONTINUATION_ATOL,
        rtol=0.0,
    )


def test_stochastic_euler_ou_stationary_mean_and_variance_match_theory() -> None:
    simulator = make_simulator(
        "ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 20.0),
        dt=0.01,
        dtmax=0.01,
        max_steps=4096,
    )
    set_repeat_ensemble_and_seed(simulator, ensemble_size=MOMENT_ENSEMBLE, seed=999)

    simulator.transient()

    final_state = simulator.get_final_state().reshape(-1)
    expected_mean = ou_stationary_mean(1.0)
    expected_variance = ou_stationary_variance(0.5)

    np.testing.assert_allclose(final_state.mean(), expected_mean, atol=0.05, rtol=0.0)
    np.testing.assert_allclose(final_state.var(), expected_variance, atol=0.02, rtol=0.0)
