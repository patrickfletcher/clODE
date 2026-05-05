from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pyopencl")

import clode
from test.core_numerics.helpers import make_simulator, set_repeat_ensemble_and_seed


FIXED_DT = 0.05
FIXED_MAX_STEPS = 256


def _make_pyopencl_simulator(monkeypatch: pytest.MonkeyPatch, **kwargs: object) -> clode.Simulator:
    monkeypatch.setenv("_CLODE_BACKEND", "pyopencl")
    return make_simulator(**kwargs)


def test_pyopencl_transient_backend_matches_cpp_for_stable_linear_rk4(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("_CLODE_BACKEND", raising=False)
    cpp = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )
    pyopencl_backend = _make_pyopencl_simulator(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    cpp.transient()
    pyopencl_backend.transient()

    np.testing.assert_allclose(
        pyopencl_backend.get_final_state(),
        cpp.get_final_state(),
        atol=1e-6,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        pyopencl_backend.get_dt(),
        cpp.get_dt(),
        atol=1e-7,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        pyopencl_backend.get_final_time(),
        cpp.get_final_time(),
        atol=1e-7,
        rtol=0.0,
    )


def test_pyopencl_transient_backend_matches_current_continuation_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    full = _make_pyopencl_simulator(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )
    split = _make_pyopencl_simulator(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 0.5),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    full.transient()

    split.transient()
    first_final_time = float(split.get_final_time()[0])
    split.set_tspan((first_final_time, 1.0))
    split.transient()

    np.testing.assert_allclose(
        split.get_final_state(),
        full.get_final_state(),
        atol=1e-6,
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
        atol=1e-7,
        rtol=0.0,
    )


def test_pyopencl_transient_backend_matches_cpp_for_seeded_stochastic_ensemble(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("_CLODE_BACKEND", raising=False)
    cpp = make_simulator(
        "ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=0.125,
        dtmax=0.125,
        max_steps=64,
    )
    pyopencl_backend = _make_pyopencl_simulator(
        monkeypatch,
        model_name="ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=0.125,
        dtmax=0.125,
        max_steps=64,
    )

    set_repeat_ensemble_and_seed(cpp, ensemble_size=64, seed=321)
    set_repeat_ensemble_and_seed(pyopencl_backend, ensemble_size=64, seed=321)

    cpp.transient()
    pyopencl_backend.transient()

    np.testing.assert_allclose(
        pyopencl_backend.get_final_state(),
        cpp.get_final_state(),
        atol=2e-7,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        pyopencl_backend.get_dt(),
        cpp.get_dt(),
        atol=2e-7,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        pyopencl_backend.get_final_time(),
        cpp.get_final_time(),
        atol=2e-7,
        rtol=0.0,
    )


def test_pyopencl_transient_backend_is_seed_reproducible_for_stochastic_ensemble(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    first = _make_pyopencl_simulator(
        monkeypatch,
        model_name="ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=0.125,
        dtmax=0.125,
        max_steps=64,
    )
    second = _make_pyopencl_simulator(
        monkeypatch,
        model_name="ornstein_uhlenbeck",
        stepper=clode.Stepper.stochastic_euler,
        t_span=(0.0, 1.0),
        dt=0.125,
        dtmax=0.125,
        max_steps=64,
    )

    set_repeat_ensemble_and_seed(first, ensemble_size=64, seed=321)
    set_repeat_ensemble_and_seed(second, ensemble_size=64, seed=321)

    first.transient()
    second.transient()

    np.testing.assert_array_equal(first.get_final_state(), second.get_final_state())
    np.testing.assert_array_equal(first.get_dt(), second.get_dt())
    np.testing.assert_array_equal(first.get_final_time(), second.get_final_time())


def test_pyopencl_transient_backend_resets_dt_after_solver_parameter_update(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    simulator = _make_pyopencl_simulator(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.1,
        dtmax=0.1,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(update_x0=False, fetch_results=False)
    simulator.set_solver_parameters(dt=0.025, dtmax=0.025)

    np.testing.assert_allclose(np.asarray(simulator.get_dt()).reshape(-1), [0.025])