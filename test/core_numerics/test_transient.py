from __future__ import annotations

import numpy as np
import pytest

import clode
from test.core_numerics.helpers import make_simulator
from test.core_numerics.reference import stable_linear_state

FIXED_DT = 0.05
FIXED_MAX_STEPS = 64
RK4_ATOL = 1e-6
DOPRI_ATOL = 2e-6


def test_rk4_stable_linear_matches_exact_final_state() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient()

    final_time = float(simulator.get_final_time()[0])
    expected = stable_linear_state(final_time)
    actual = simulator.get_final_state()[0]

    np.testing.assert_allclose(actual, expected, atol=RK4_ATOL, rtol=0.0)


def test_dormand_prince_stable_linear_matches_exact_final_state() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.dormand_prince,
        t_span=(0.0, 1.0),
        dt=0.05,
        dtmax=0.1,
        abstol=1e-7,
        reltol=1e-6,
        max_steps=256,
    )

    simulator.transient()

    final_time = float(simulator.get_final_time()[0])
    expected = stable_linear_state(final_time)
    actual = simulator.get_final_state()[0]

    assert final_time == pytest.approx(1.0, abs=5e-6)
    np.testing.assert_allclose(actual, expected, atol=DOPRI_ATOL, rtol=0.0)


def test_rk4_large_origin_fixed_step_reports_counter_reconstructed_final_time() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(1_000_000.0, 1_000_100.0),
        dt=0.01,
        dtmax=0.01,
        max_steps=20_000,
    )

    simulator.transient(update_x0=False, fetch_results=False)

    final_time = float(simulator.get_final_time()[0])

    assert final_time > 1_000_050.0
    assert final_time == pytest.approx(float(np.float32(1_000_100.0)), abs=0.0)


def test_rk4_stable_linear_ensemble_matches_exact_final_state() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.set_ensemble(
        variables={
            "x": np.array([2.0, 1.0, -1.0, 0.5], dtype=np.float64),
            "y": np.array([-1.5, 0.25, -0.5, 3.0], dtype=np.float64),
        },
        parameters={
            "a": np.array([0.5, 0.25, 1.0, 0.75], dtype=np.float64),
            "b": np.array([1.25, 0.5, 0.75, 1.5], dtype=np.float64),
        },
    )

    simulator.transient()

    final_time = float(simulator.get_final_time().reshape(-1)[0])
    variables = np.array(
        [[2.0, -1.5], [1.0, 0.25], [-1.0, -0.5], [0.5, 3.0]], dtype=np.float64
    )
    parameters = np.array(
        [[0.5, 1.25], [0.25, 0.5], [1.0, 0.75], [0.75, 1.5]], dtype=np.float64
    )
    expected = np.column_stack(
        (
            variables[:, 0] * np.exp(-parameters[:, 0] * final_time),
            variables[:, 1] * np.exp(-parameters[:, 1] * final_time),
        )
    )

    np.testing.assert_allclose(
        simulator.get_final_state(),
        expected,
        atol=RK4_ATOL,
        rtol=0.0,
    )


@pytest.mark.parametrize(
    ("stepper", "dt", "dtmax", "split_end", "full_end", "kwargs", "atol"),
    [
        (clode.Stepper.rk4, FIXED_DT, FIXED_DT, 0.5, 1.0, {}, RK4_ATOL),
        (
            clode.Stepper.dormand_prince,
            0.05,
            0.1,
            0.6,
            1.0,
            {"abstol": 1e-7, "reltol": 1e-6},
            DOPRI_ATOL,
        ),
    ],
)
def test_deterministic_transient_continuation_matches_single_run(
    stepper: clode.Stepper,
    dt: float,
    dtmax: float,
    split_end: float,
    full_end: float,
    kwargs: dict[str, float],
    atol: float,
) -> None:
    full = make_simulator(
        "stable_linear",
        stepper=stepper,
        t_span=(0.0, full_end),
        dt=dt,
        dtmax=dtmax,
        max_steps=256,
        **kwargs,
    )
    full.transient()

    split = make_simulator(
        "stable_linear",
        stepper=stepper,
        t_span=(0.0, split_end),
        dt=dt,
        dtmax=dtmax,
        max_steps=256,
        **kwargs,
    )
    split.transient()
    first_final_time = float(split.get_final_time()[0])
    split.set_tspan((first_final_time, full_end))
    split.transient()

    np.testing.assert_allclose(
        split.get_final_state(),
        full.get_final_state(),
        atol=atol,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        split.get_dt(),
        full.get_dt(),
        atol=1e-7,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        split.get_final_time(),
        full.get_final_time(),
        atol=0.0,
        rtol=0.0,
    )
