from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pyopencl")

import clode
from test.core_numerics.helpers import make_trajectory_simulator, structured_to_array


FIXED_DT = 0.05
FIXED_MAX_STEPS = 64


def _make_pyopencl_trajectory(
    monkeypatch: pytest.MonkeyPatch, **kwargs: object
) -> clode.TrajectorySimulator:
    monkeypatch.setenv("_CLODE_BACKEND", "pyopencl")
    return make_trajectory_simulator(**kwargs)


def test_pyopencl_trajectory_backend_matches_cpp_for_stable_linear_samples(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("_CLODE_BACKEND", raising=False)
    cpp = make_trajectory_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=FIXED_MAX_STEPS,
        nout=1,
    )
    pyopencl_backend = _make_pyopencl_trajectory(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=FIXED_MAX_STEPS,
        nout=1,
    )

    cpp_result = cpp.trajectory()
    pyopencl_result = pyopencl_backend.trajectory()

    np.testing.assert_allclose(
        np.asarray(pyopencl_result.t, dtype=np.float64),
        np.asarray(cpp_result.t, dtype=np.float64),
        atol=1e-7,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        structured_to_array(pyopencl_result, "x"),
        structured_to_array(cpp_result, "x"),
        atol=1e-6,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        structured_to_array(pyopencl_result, "dx"),
        structured_to_array(cpp_result, "dx"),
        atol=1e-6,
        rtol=0.0,
    )


def test_pyopencl_trajectory_backend_matches_cpp_for_aux_samples(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("_CLODE_BACKEND", raising=False)
    cpp = make_trajectory_simulator(
        "stable_linear_aux",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=FIXED_MAX_STEPS,
        nout=1,
    )
    pyopencl_backend = _make_pyopencl_trajectory(
        monkeypatch,
        model_name="stable_linear_aux",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=FIXED_MAX_STEPS,
        nout=1,
    )

    cpp_result = cpp.trajectory()
    pyopencl_result = pyopencl_backend.trajectory()

    np.testing.assert_allclose(
        structured_to_array(pyopencl_result, "aux"),
        structured_to_array(cpp_result, "aux"),
        atol=1e-6,
        rtol=0.0,
    )


def test_pyopencl_trajectory_backend_respects_nout_and_max_store_contract(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    simulator = _make_pyopencl_trajectory(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.25,
        dtmax=0.25,
        max_steps=32,
        max_store=3,
        nout=2,
    )

    result = simulator.trajectory()

    np.testing.assert_allclose(
        np.asarray(result.t, dtype=np.float64),
        np.asarray([0.0, 0.5, 1.0], dtype=np.float64),
        atol=0.0,
        rtol=0.0,
    )
    assert int(simulator._device_n_stored[0]) == 2


def test_pyopencl_trajectory_backend_reallocates_when_max_store_changes(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    simulator = _make_pyopencl_trajectory(
        monkeypatch,
        model_name="stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 0.5),
        dt=0.25,
        dtmax=0.25,
        max_steps=32,
        max_store=3,
        nout=1,
    )

    first = simulator.trajectory(update_x0=False)
    assert len(np.asarray(first.t, dtype=np.float64)) == 3
    assert len(simulator._integrator.get_t()) == 3

    simulator.set_solver_parameters(max_store=5)
    second = simulator.trajectory(update_x0=False)

    assert len(np.asarray(second.t, dtype=np.float64)) == 3
    assert len(simulator._integrator.get_t()) == 5
    assert int(simulator._device_n_stored[0]) == 2