import pytest
import numpy as np

import clode
from test.core_numerics.helpers import (
    HOPF_AUX,
    HOPF_PARAMETERS,
    HOPF_VARIABLES,
    STABLE_LINEAR_AUX,
    STABLE_LINEAR_PARAMETERS,
    STABLE_LINEAR_VARIABLES,
    device_kwargs_for_tests,
    make_feature_simulator,
    make_simulator,
    model_path,
)
from test.core_numerics.reference import fixed_step_step_count


FIXED_DT = 0.05
FIXED_MAX_STEPS = 512


def test_get_tspan_returns_current_device_window() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(1.25, 2.5),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    assert simulator.get_tspan() == (1.25, 2.5)


def test_get_final_state_can_be_fetched_twice() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(fetch_results=False)
    first = simulator.get_final_state().copy()
    second = simulator.get_final_state()

    np.testing.assert_allclose(second, first)


def test_set_solver_parameters_resets_device_dt() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.1,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(update_x0=False, fetch_results=False)
    simulator.set_solver_parameters(dt=0.025, dtmax=0.025)

    np.testing.assert_allclose(np.asarray(simulator.get_dt()).reshape(-1), [0.025])


def test_get_initial_state_after_update_x0_pulls_runtime_state_back_into_ivp() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(update_x0=True, fetch_results=False)

    final_state = simulator.get_final_state()
    next_initial_state = simulator.get_initial_state()

    np.testing.assert_allclose(next_initial_state, final_state)
    np.testing.assert_allclose(simulator.ivp.get_initial_state(), final_state)


def test_set_tspan_invalidates_previous_transient_results() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(update_x0=False, fetch_results=False)
    simulator.get_final_state()
    simulator.get_final_time()

    simulator.set_tspan((1.0, 2.0))

    with pytest.raises(ValueError, match="final state"):
        simulator.get_final_state()

    with pytest.raises(ValueError, match="final time"):
        simulator.get_final_time()


def test_set_tspan_invalidates_cached_feature_results() -> None:
    simulator = make_feature_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    first = simulator.features(update_x0=False)

    assert first is not None

    simulator.set_tspan((0.0, 2.0))

    with pytest.raises(ValueError, match=r"features\(\)"):
        simulator.get_observer_results()


def test_feature_observer_switch_rebuilds_and_runs() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.basic_all_variables,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    initial_feature_names = simulator.get_feature_names()

    simulator.set_observer(clode.Observer.basic)
    result = simulator.features(update_x0=False)

    assert result is not None
    assert len(initial_feature_names) > len(result.get_feature_names())
    assert result.get_feature_names() == simulator.get_feature_names()
    assert result.get_feature_names() == [
        "max x",
        "min x",
        "mean x",
        "max dx/dt",
        "min dx/dt",
        "step count",
    ]


def test_features_initialize_observer_path_refreshes_results() -> None:
    simulator = make_feature_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    first = simulator.features(update_x0=False)
    second = simulator.features(
        t_span=(0.0, 2.0),
        initialize_observer=True,
        update_x0=False,
    )

    assert first is not None
    assert second is not None
    assert int(first.get_var_count("step")) == fixed_step_step_count(0.0, 1.0, FIXED_DT)
    assert int(second.get_var_count("step")) == fixed_step_step_count(0.0, 2.0, FIXED_DT)


def test_feature_max_event_timestamps_rebuild_updates_event_storage() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.local_max,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 8.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    simulator.set_observer_parameters(max_event_timestamps=3)
    result = simulator.features(update_x0=False)

    assert result is not None
    feature_names = result.get_feature_names()
    assert len([name for name in feature_names if name.startswith("localmax event time")]) == 3
    assert "-DN_STORE_EVENTS=3" in simulator._integrator.get_program_string()

