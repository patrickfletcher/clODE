import inspect

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
    make_trajectory_simulator,
    model_path,
)
from test.core_numerics.reference import fixed_step_step_count, fixed_step_time_grid


FIXED_DT = 0.05
FIXED_MAX_STEPS = 512
IN_LOOP_NO_PROGRESS_TSPAN = (1.0e6, 1.0e6 + 0.2)
IN_LOOP_NO_PROGRESS_DT = 0.01


@pytest.mark.parametrize(
    "constructor",
    [clode.Simulator, clode.TrajectorySimulator, clode.FeatureSimulator],
)
def test_simulator_constructor_defaults_match_solver_params_defaults(
    constructor: type[clode.Simulator],
) -> None:
    signature = inspect.signature(constructor.__init__)
    defaults = clode.SolverParams()

    assert signature.parameters["dt"].default == defaults.dt
    assert signature.parameters["dtmax"].default == defaults.dtmax
    assert signature.parameters["abstol"].default == defaults.abstol
    assert signature.parameters["reltol"].default == defaults.reltol
    assert signature.parameters["max_steps"].default == defaults.max_steps
    assert signature.parameters["max_store"].default == defaults.max_store
    assert signature.parameters["nout"].default == defaults.nout


def test_scalar_solver_args_match_explicit_solver_params_bundle() -> None:
    bundle = clode.SolverParams(
        dt=0.125,
        dtmax=0.25,
        abstol=1e-7,
        reltol=1e-5,
        max_steps=2048,
        max_store=64,
        nout=4,
    )

    explicit = clode.Simulator(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        solver_parameters=bundle,
        single_precision=True,
        **device_kwargs_for_tests(),
    )
    scalar = clode.Simulator(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        dt=bundle.dt,
        dtmax=bundle.dtmax,
        abstol=bundle.abstol,
        reltol=bundle.reltol,
        max_steps=bundle.max_steps,
        max_store=bundle.max_store,
        nout=bundle.nout,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert explicit.get_solver_parameters() == bundle
    assert scalar.get_solver_parameters() == bundle


def test_simulator_copies_solver_parameter_bundles() -> None:
    initial = clode.SolverParams(dt=0.2, dtmax=0.3, max_steps=128)
    replacement = clode.SolverParams(dt=0.05, dtmax=0.06, max_steps=64)

    simulator = clode.Simulator(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        solver_parameters=initial,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    simulator.set_solver_parameters(dt=0.125)
    assert initial.dt == 0.2

    simulator.set_solver_parameters(solver_parameters=replacement)
    simulator.set_solver_parameters(dt=0.03125)

    assert replacement.dt == 0.05
    assert simulator.get_solver_parameters().dt == pytest.approx(0.03125)


def test_get_tspan_returns_current_device_window() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(1.25, 2.5),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    assert simulator.get_tspan() == (1.25, 2.5)


def test_advance_tspan_to_attained_final_time_uses_attained_fixed_step_end() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.3,
        dtmax=0.3,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(update_x0=False, fetch_results=False)
    attained_final_time = float(simulator.get_final_time().reshape(-1)[0])
    expected_attained = float(fixed_step_time_grid(0.0, 1.0, 0.3)[-1])

    assert attained_final_time == expected_attained
    assert simulator.advance_tspan_to_attained_final_time() == (expected_attained, expected_attained + 1.0)
    assert simulator.get_tspan() == (expected_attained, expected_attained + 1.0)


def test_shift_tspan_keeps_requested_window_semantics_after_fixed_step_run() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.3,
        dtmax=0.3,
        max_steps=FIXED_MAX_STEPS,
    )

    simulator.transient(update_x0=False, fetch_results=False)
    simulator.shift_tspan()

    assert simulator.get_tspan() == pytest.approx((1.0, 2.0))


def test_advance_tspan_to_attained_final_time_rejects_diverged_final_times() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )
    simulator.set_repeat_ensemble(2)
    simulator._solver_state.final_time = np.array([0.5, 0.75], dtype=np.float64)

    with pytest.raises(ValueError, match="different final times"):
        simulator.advance_tspan_to_attained_final_time()


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


def test_get_status_reports_max_steps_exhaustion() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.3,
        dtmax=0.3,
        max_steps=1,
    )

    simulator.transient(update_x0=False, fetch_results=False)

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.MAX_STEPS_REACHED],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [1])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.3])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [0.3])


def test_get_status_reports_no_progress_when_float32_time_cannot_advance() -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(1.0e6, 1.0e6 + 0.01),
        dt=0.001,
        dtmax=0.001,
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    simulator.transient(update_x0=False, fetch_results=False)

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.NO_PROGRESS],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [0])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.0])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [1.0e6])


def test_trajectory_status_reports_no_progress_when_float32_time_cannot_advance() -> None:
    simulator = make_trajectory_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(1.0e6, 1.0e6 + 0.01),
        dt=0.001,
        dtmax=0.001,
        max_steps=FIXED_MAX_STEPS,
        max_store=32,
        nout=1,
        single_precision=True,
    )

    trajectory = simulator.trajectory(update_x0=False)
    output = trajectory if not isinstance(trajectory, list) else trajectory[0]

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.NO_PROGRESS],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [0])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.0])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [1.0e6])
    assert len(output.t) == 1


def test_feature_status_reports_no_progress_when_float32_time_cannot_advance() -> None:
    simulator = make_feature_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(1.0e6, 1.0e6 + 0.01),
        dt=0.001,
        dtmax=0.001,
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    simulator.features(update_x0=False)

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.NO_PROGRESS],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [0])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.0])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [1.0e6])


@pytest.mark.parametrize(
    "stepper",
    [clode.Stepper.rk4, clode.Stepper.dormand_prince],
)
def test_get_status_reports_no_progress_for_in_loop_time_stall(
    stepper: clode.Stepper,
) -> None:
    simulator = make_simulator(
        "stable_linear",
        stepper=stepper,
        t_span=IN_LOOP_NO_PROGRESS_TSPAN,
        dt=IN_LOOP_NO_PROGRESS_DT,
        dtmax=IN_LOOP_NO_PROGRESS_DT,
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    simulator.transient(update_x0=False, fetch_results=False)

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.NO_PROGRESS],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [0])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.0])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [1.0e6])


@pytest.mark.parametrize(
    "stepper",
    [clode.Stepper.rk4, clode.Stepper.dormand_prince],
)
def test_trajectory_status_reports_no_progress_for_in_loop_time_stall(
    stepper: clode.Stepper,
) -> None:
    simulator = make_trajectory_simulator(
        "stable_linear",
        stepper=stepper,
        t_span=IN_LOOP_NO_PROGRESS_TSPAN,
        dt=IN_LOOP_NO_PROGRESS_DT,
        dtmax=IN_LOOP_NO_PROGRESS_DT,
        max_steps=FIXED_MAX_STEPS,
        max_store=32,
        nout=1,
        single_precision=True,
    )

    trajectory = simulator.trajectory(update_x0=False)
    output = trajectory if not isinstance(trajectory, list) else trajectory[0]

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.NO_PROGRESS],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [0])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.0])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [1.0e6])
    assert len(output.t) == 1


@pytest.mark.parametrize(
    "stepper",
    [clode.Stepper.rk4, clode.Stepper.dormand_prince],
)
def test_feature_status_reports_no_progress_for_in_loop_time_stall(
    stepper: clode.Stepper,
) -> None:
    simulator = make_feature_simulator(
        "stable_linear",
        stepper=stepper,
        t_span=IN_LOOP_NO_PROGRESS_TSPAN,
        dt=IN_LOOP_NO_PROGRESS_DT,
        dtmax=IN_LOOP_NO_PROGRESS_DT,
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    simulator.features(update_x0=False)

    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.NO_PROGRESS],
    )
    np.testing.assert_array_equal(simulator.get_step_count().reshape(-1), [0])
    np.testing.assert_allclose(simulator.get_last_accepted_dt().reshape(-1), [0.0])
    np.testing.assert_allclose(simulator.get_final_time().reshape(-1), [1.0e6])


def test_feature_status_reports_terminal_event_stop() -> None:
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
    simulator.set_observer_parameters(max_event_count=1, max_event_timestamps=1)

    result = simulator.features(update_x0=False)

    assert result is not None
    assert int(result.get_var_count("event")) == 1
    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.TERMINAL_EVENT_REACHED],
    )
    assert float(simulator.get_final_time().reshape(-1)[0]) < 8.0


def test_feature_status_reports_terminal_event_stop_for_neighborhood_return() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.normalized_neighborhood_return,
        observer_configuration=clode.NeighborhoodReturnConfig(
            event_var="x",
            feature_var="x",
            anchor_threshold=0.35,
            radius=0.2,
            max_event_count=1,
            max_event_timestamps=1,
        ),
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 8.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    result = simulator.features(update_x0=False)

    assert result is not None
    assert int(result.get_var_count("event")) == 1
    np.testing.assert_array_equal(
        simulator.get_status().reshape(-1),
        [clode.SolverStatus.TERMINAL_EVENT_REACHED],
    )
    assert float(simulator.get_final_time().reshape(-1)[0]) < 8.0


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


def test_output_only_trajectory_config_change_preserves_solver_state() -> None:
    simulator = make_trajectory_simulator(
        "stable_linear",
        stepper=clode.Stepper.dormand_prince,
        t_span=(0.0, 4.0),
        dt=0.01,
        dtmax=0.2,
        abstol=1e-9,
        reltol=1e-8,
        max_steps=1024,
        max_store=16,
        nout=1,
    )

    simulator.trajectory(update_x0=False, fetch_results=False)

    current_dt = simulator.get_dt().copy()
    final_time = simulator.get_final_time().copy()

    assert current_dt.reshape(-1)[0] > 0.02

    simulator.set_solver_parameters(max_store=4, nout=2)

    np.testing.assert_allclose(simulator.get_dt(), current_dt)
    np.testing.assert_allclose(simulator.get_final_time(), final_time)

    with pytest.raises(ValueError, match="trajectory data"):
        simulator.get_trajectory()


def test_feature_observer_switch_rebuilds_and_runs() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.summary,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    initial_feature_names = simulator.get_feature_names()

    simulator.set_observer(clode.Observer.local_max)
    result = simulator.features(update_x0=False)

    assert result is not None
    assert result.get_feature_names() == simulator.get_feature_names()
    assert initial_feature_names != result.get_feature_names()


def test_feature_simulator_constructor_defaults_match_observer_params() -> None:
    simulator = make_feature_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    assert simulator.get_observer_parameters() == clode.ObserverParams()


def test_feature_simulator_threshold_crossing_configuration_round_trips() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.summary,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )
    config = clode.ThresholdCrossingConfig(
        event_var="y",
        feature_var="x",
        threshold=0.75,
        direction=clode.EventDirection.falling,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=2,
    )

    simulator.set_observer_configuration(
        config,
        observer=clode.Observer.threshold_crossing,
    )

    assert simulator.get_observer_configuration() == config
    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=7,
        event_direction=clode.EventDirection.falling,
        max_event_timestamps=2,
        min_amp=0.125,
        x_up_threshold=0.75,
        x_down_threshold=0.75,
    )


def test_feature_simulator_normalized_threshold_crossing_configuration_round_trips() -> None:
    config = clode.ThresholdCrossingConfig(
        event_var="y",
        feature_var="x",
        threshold=0.6,
        direction=clode.EventDirection.falling,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
    )
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.normalized_threshold_crossing,
        observer_configuration=config,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator.get_observer_configuration() == config
    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=7,
        event_direction=clode.EventDirection.falling,
        max_event_timestamps=4,
        min_amp=0.125,
        x_up_threshold=0.6,
        x_down_threshold=0.6,
    )


def test_feature_simulator_normalized_schmitt_trigger_configuration_round_trips() -> None:
    config = clode.SchmittTriggerConfig(
        event_var="y",
        feature_var="x",
        x_up_threshold=0.6,
        x_down_threshold=0.4,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
    )
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.normalized_schmitt_trigger,
        observer_configuration=config,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator.get_observer_configuration() == config
    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=7,
        event_direction=clode.EventDirection.either,
        max_event_timestamps=4,
        min_amp=0.125,
        x_up_threshold=0.6,
        x_down_threshold=0.4,
    )


def test_feature_simulator_absolute_schmitt_trigger_configuration_round_trips() -> None:
    config = clode.SchmittTriggerConfig(
        event_var="y",
        feature_var="x",
        x_up_threshold=0.6,
        x_down_threshold=0.4,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
    )
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.schmitt_trigger,
        observer_configuration=config,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator.get_observer_configuration() == config
    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=7,
        event_direction=clode.EventDirection.either,
        max_event_timestamps=4,
        min_amp=0.125,
        x_up_threshold=0.6,
        x_down_threshold=0.4,
    )


def test_schmitt_trigger_config_allows_equal_thresholds() -> None:
    config = clode.SchmittTriggerConfig(
        event_var="x",
        x_up_threshold=0.4,
        x_down_threshold=0.4,
    )

    assert config.x_up_threshold == pytest.approx(0.4)
    assert config.x_down_threshold == pytest.approx(0.4)


def test_schmitt_trigger_config_rejects_inverted_thresholds() -> None:
    with pytest.raises(ValueError, match="x_up_threshold >= x_down_threshold"):
        clode.SchmittTriggerConfig(
            event_var="x",
            x_up_threshold=0.3,
            x_down_threshold=0.4,
        )


def test_feature_simulator_rejects_dx_thresholds_for_semantic_schmitt() -> None:
    with pytest.raises(ValueError, match="do not support dx thresholds"):
        clode.FeatureSimulator(
            src_file=model_path("hopf_normal_form.cl"),
            variables=HOPF_VARIABLES.copy(),
            parameters=HOPF_PARAMETERS.copy(),
            aux=HOPF_AUX.copy(),
            num_noise=0,
            observer=clode.Observer.schmitt_trigger,
            event_var="y",
            feature_var="x",
            observer_x_up_thresh=0.6,
            observer_x_down_thresh=0.4,
            observer_dx_up_thresh=0.2,
            stepper=clode.Stepper.rk4,
            dt=FIXED_DT,
            dtmax=FIXED_DT,
            t_span=(0.0, 1.0),
            max_steps=FIXED_MAX_STEPS,
            single_precision=True,
            **device_kwargs_for_tests(),
        )


def test_feature_simulator_rejects_out_of_range_normalized_thresholds() -> None:
    with pytest.raises(ValueError, match="x_up_threshold must be between 0.0 and 1.0 inclusive"):
        clode.FeatureSimulator(
            src_file=model_path("hopf_normal_form.cl"),
            variables=HOPF_VARIABLES.copy(),
            parameters=HOPF_PARAMETERS.copy(),
            aux=HOPF_AUX.copy(),
            num_noise=0,
            observer=clode.Observer.normalized_schmitt_trigger,
            observer_configuration=clode.SchmittTriggerConfig(
                event_var="y",
                feature_var="x",
                x_up_threshold=1.1,
                x_down_threshold=0.4,
            ),
            stepper=clode.Stepper.rk4,
            dt=FIXED_DT,
            dtmax=FIXED_DT,
            t_span=(0.0, 1.0),
            max_steps=FIXED_MAX_STEPS,
            single_precision=True,
            **device_kwargs_for_tests(),
        )

    with pytest.raises(ValueError, match="threshold must be between 0.0 and 1.0 inclusive"):
        clode.FeatureSimulator(
            src_file=model_path("hopf_normal_form.cl"),
            variables=HOPF_VARIABLES.copy(),
            parameters=HOPF_PARAMETERS.copy(),
            aux=HOPF_AUX.copy(),
            num_noise=0,
            observer=clode.Observer.normalized_threshold_crossing,
            observer_configuration=clode.ThresholdCrossingConfig(
                event_var="y",
                threshold=1.1,
            ),
            stepper=clode.Stepper.rk4,
            dt=FIXED_DT,
            dtmax=FIXED_DT,
            t_span=(0.0, 1.0),
            max_steps=FIXED_MAX_STEPS,
            single_precision=True,
            **device_kwargs_for_tests(),
        )


def test_neighborhood_return_config_rejects_invalid_normalized_inputs() -> None:
    with pytest.raises(ValueError, match="anchor_threshold must be between 0.0 and 1.0 inclusive"):
        clode.NeighborhoodReturnConfig(anchor_threshold=1.1)

    with pytest.raises(ValueError, match="radius must be greater than 0.0"):
        clode.NeighborhoodReturnConfig(radius=0.0)


def test_feature_simulator_local_maximum_configuration_round_trips() -> None:
    config = clode.LocalMaximumConfig(
        event_var="y",
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
    )
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer_configuration=config,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator.get_observer_configuration() == config
    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=1,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
    )


def test_feature_simulator_neighborhood_return_configuration_round_trips() -> None:
    config = clode.NeighborhoodReturnConfig(
        event_var="y",
        feature_var="x",
        anchor_threshold=0.35,
        radius=0.2,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
    )
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer_configuration=config,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator.get_observer_configuration() == config
    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        min_amp=0.125,
        max_event_count=7,
        max_event_timestamps=4,
        nhood_radius=0.2,
        x_down_threshold=0.35,
    )


def test_feature_simulator_requires_explicit_observer_for_shared_observer_configuration() -> None:
    with pytest.raises(ValueError, match="requires an explicit observer"):
        clode.FeatureSimulator(
            src_file=model_path("hopf_normal_form.cl"),
            variables=HOPF_VARIABLES.copy(),
            parameters=HOPF_PARAMETERS.copy(),
            aux=HOPF_AUX.copy(),
            num_noise=0,
            observer_configuration=clode.ThresholdCrossingConfig(event_var="x"),
            stepper=clode.Stepper.rk4,
            dt=FIXED_DT,
            dtmax=FIXED_DT,
            t_span=(0.0, 1.0),
            max_steps=FIXED_MAX_STEPS,
            single_precision=True,
            **device_kwargs_for_tests(),
        )


def test_feature_simulator_rejects_mixed_observer_configuration_and_compatibility_args() -> None:
    with pytest.raises(ValueError, match="observer_configuration cannot be combined"):
        clode.FeatureSimulator(
            src_file=model_path("hopf_normal_form.cl"),
            variables=HOPF_VARIABLES.copy(),
            parameters=HOPF_PARAMETERS.copy(),
            aux=HOPF_AUX.copy(),
            num_noise=0,
            observer_configuration=clode.ThresholdCrossingConfig(event_var="x"),
            event_var="x",
            stepper=clode.Stepper.rk4,
            dt=FIXED_DT,
            dtmax=FIXED_DT,
            t_span=(0.0, 1.0),
            max_steps=FIXED_MAX_STEPS,
            single_precision=True,
            **device_kwargs_for_tests(),
        )


def test_feature_simulator_copies_observer_parameter_bundles() -> None:
    initial = clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=12,
        max_event_timestamps=3,
        min_amp=0.25,
    )
    replacement = clode.ObserverParams(
        e_var_ix=0,
        f_var_ix=1,
        max_event_count=7,
        max_event_timestamps=1,
        min_amp=0.5,
    )

    simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.summary,
        observer_parameters=initial,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    initial.min_amp = 99.0
    initial.max_event_timestamps = 99

    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=12,
        max_event_timestamps=3,
        min_amp=0.25,
    )

    simulator.set_observer_parameters(op=replacement)
    replacement.f_var_ix = 0
    replacement.max_event_count = 99

    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=0,
        f_var_ix=1,
        max_event_count=7,
        max_event_timestamps=1,
        min_amp=0.5,
    )


def test_feature_simulator_resolves_scalar_observer_args_from_variable_names() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.schmitt_trigger,
        event_var="y",
        feature_var="x",
        observer_max_event_count=7,
        observer_max_event_timestamps=4,
        observer_min_x_amp=0.125,
        observer_min_imi=0.25,
        observer_neighbourhood_radius=0.5,
        observer_x_up_thresh=0.6,
        observer_x_down_thresh=0.4,
        observer_eps_dx=0.05,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator.get_observer_parameters() == clode.ObserverParams(
        e_var_ix=1,
        f_var_ix=0,
        max_event_count=7,
        max_event_timestamps=4,
        min_amp=0.125,
        min_imi=0.25,
        nhood_radius=0.5,
        x_up_threshold=0.6,
        x_down_threshold=0.4,
        eps_dx=0.05,
    )


def test_feature_simulator_rejects_unknown_observer_variable_names() -> None:
    with pytest.raises(ValueError, match="Unknown event_var 'z'"):
        clode.FeatureSimulator(
            src_file=model_path("stable_linear.cl"),
            variables=STABLE_LINEAR_VARIABLES.copy(),
            parameters=STABLE_LINEAR_PARAMETERS.copy(),
            observer=clode.Observer.summary,
            event_var="z",
            stepper=clode.Stepper.rk4,
            dt=FIXED_DT,
            dtmax=FIXED_DT,
            t_span=(0.0, 1.0),
            max_steps=FIXED_MAX_STEPS,
            single_precision=True,
            **device_kwargs_for_tests(),
        )


def test_features_initialize_observer_path_refreshes_results() -> None:
    simulator = make_feature_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    first = simulator.features(update_x0=False)
    first_step_count = int(simulator.get_step_count().reshape(-1)[0])
    second = simulator.features(
        t_span=(0.0, 2.0),
        initialize_observer=True,
        update_x0=False,
    )
    second_step_count = int(simulator.get_step_count().reshape(-1)[0])

    assert first is not None
    assert second is not None
    assert first_step_count == fixed_step_step_count(0.0, 1.0, FIXED_DT)
    assert second_step_count == fixed_step_step_count(0.0, 2.0, FIXED_DT)


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
    assert len([name for name in feature_names if name.startswith("localmax time")]) == 3
    assert "-DN_STORE_EVENTS=3" in simulator._integrator.get_program_string()


def test_feature_max_event_timestamps_rebuild_updates_local_maximum_storage() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.local_max,
        observer_configuration=clode.LocalMaximumConfig(
            event_var="x",
            max_event_count=8,
            max_event_timestamps=3,
        ),
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 8.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    result = simulator.features(update_x0=False)

    assert result is not None
    feature_names = result.get_feature_names()
    assert len(
        [
            name
            for name in feature_names
            if name.startswith("localmax time")
        ]
    ) == 3
    assert "-DN_STORE_EVENTS=3" in simulator._integrator.get_program_string()


def test_feature_var_change_rebuilds_program_but_reuses_feature_buffers() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.summary,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    first = simulator.features(update_x0=False)
    program_bundle = simulator._integrator._program_bundle
    feature_buffers = simulator._integrator._feature_buffers

    assert first is not None
    assert program_bundle is not None
    assert feature_buffers is not None
    assert first.get_feature_names()[0] == "max x"

    simulator.set_observer_parameters(feature_var="y")

    assert simulator._integrator._program_bundle is program_bundle
    assert simulator._integrator._feature_buffers is feature_buffers

    with pytest.raises(ValueError, match=r"features\(\)"):
        simulator.get_observer_results()

    second = simulator.features(update_x0=False)

    assert second is not None
    assert second.get_feature_names()[0] == "max x"


def test_feature_summary_selection_matches_basicall_subset() -> None:
    selection = clode.SummaryObserverSelection(
        state={"y": ("max", "mean")},
        slope={"x": "min"},
        aux={STABLE_LINEAR_AUX[0]: "mean"},
    )

    all_simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.summary,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )
    summary_simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.summary,
        summary_selection=selection,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    all_result = all_simulator.features(update_x0=False)
    summary_result = summary_simulator.features(update_x0=False)

    assert all_result is not None
    assert summary_result is not None
    assert summary_simulator.get_summary_selection() == selection
    assert summary_result.get_feature_names() == [
        "max y",
        "mean y",
        "min dx/dt",
        f"mean {STABLE_LINEAR_AUX[0]}",
    ]
    np.testing.assert_allclose(
        summary_result.get_var_max("y"),
        all_result.get_var_max("y"),
    )
    np.testing.assert_allclose(
        summary_result.get_var_mean("y"),
        all_result.get_var_mean("y"),
    )
    np.testing.assert_allclose(
        summary_result.get_var_min_slope("x"),
        all_result.get_var_min_slope("x"),
    )
    np.testing.assert_allclose(
        summary_result.get_var_mean(STABLE_LINEAR_AUX[0]),
        all_result.get_var_mean(STABLE_LINEAR_AUX[0]),
    )


def test_feature_observer_parameter_update_rejects_unknown_variable_names() -> None:
    simulator = make_feature_simulator(
        "stable_linear_aux",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=FIXED_DT,
        max_steps=FIXED_MAX_STEPS,
    )

    with pytest.raises(ValueError, match="Unknown feature_var 'z'"):
        simulator.set_observer_parameters(feature_var="z")


def test_threshold_crossing_always_exposes_event_times_and_count() -> None:
    simulator = clode.FeatureSimulator(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        observer=clode.Observer.threshold_crossing,
        observer_configuration=clode.ThresholdCrossingConfig(
            event_var="x",
            threshold=0.5,
            direction=clode.EventDirection.either,
            min_amp=0.5,
            max_event_count=8,
            max_event_timestamps=3,
        ),
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        **device_kwargs_for_tests(),
    )

    initial = simulator.features(update_x0=False)

    assert initial is not None
    feature_names = initial.get_feature_names()
    assert feature_names[:10] == [
        "event time 0",
        "event x 0",
        "event y 0",
        "event time 1",
        "event x 1",
        "event y 1",
        "event time 2",
        "event x 2",
        "event y 2",
        "event count",
    ]
    assert "period max" in feature_names
    assert "amplitude mean" in feature_names
    assert int(initial.get_var_count("event")) == 0

