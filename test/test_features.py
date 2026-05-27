from math import cos, exp, pi, sqrt
from typing import List

import numpy as np
import pytest

import clode
from test.core_numerics.helpers import device_kwargs_for_tests


def sine_curve(
    t: float,
    x_: List[float],
    p_: List[float],
    dx_: List[float],
    aux_: List[float],
    w_: List[float],
) -> None:
    x: float = x_[0]
    dilation: float = p_[0]
    dx: float = cos(t * dilation)
    dx_[0] = dx


def test_sine_curve_timestamps():
    "Test that the active timestamps of a sine curve are correct"

    # Define the parameters
    parameters = {
        "dilation": 1,
    }

    # Define the initial state
    variables = {
        "x": 0,
    }

    # Activate at t=pi/4, deactivate on the descending x_down threshold crossing.
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables=variables,
        parameters=parameters,
        aux=["dx"],
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.rk4,
        dtmax=0.001,
        dt=0.001,
        t_span=(0.0, 4 * pi),
        event_var="x",
        feature_var="x",
        observer_min_x_amp=0.5,
        observer_x_up_thresh=(2 + sqrt(2)) / 4,
        observer_x_down_thresh=0.001,
        observer_dx_down_thresh=0.001,
        observer_dx_up_thresh=0.001,
        observer_max_event_count=100,
        observer_max_event_timestamps=3,
        **device_kwargs_for_tests(),
    )

    # Run the simulation
    output = feature_simulator.features()

    assert output is not None
    event_count = int(output.get_var_count("event"))

    up_times = output.get_event_data("up")
    down_times = output.get_event_data("down")

    assert len(up_times) == 2
    assert len(down_times) == 2
    assert event_count == 2

    assert np.isclose(up_times[0], pi / 4, atol=2e-3)
    assert np.isclose(up_times[1], 9 * pi / 4, atol=2e-3)

    x_down = -1.0 + 2.0e-3
    assert np.isclose(down_times[0], pi - np.arcsin(x_down), atol=2e-3)
    assert np.isclose(down_times[1], 3 * pi - np.arcsin(x_down), atol=2e-3)

    up_times_timestamps = output.get_timestamps("up")
    down_times_timestamps = output.get_timestamps("down")

    assert len(up_times_timestamps) == len(up_times)
    assert len(down_times_timestamps) == len(down_times)

    for index in range(len(up_times)):
        assert up_times_timestamps[index] == up_times[index]

    for index in range(len(down_times)):
        assert down_times_timestamps[index] == down_times[index]


def test_threshold_2_coarse_timestamps_use_inverse_linear_interpolation() -> None:
    coarse_dt = pi / 6
    coarse_x_down_thresh = 0.25

    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0},
        parameters={"dilation": 1},
        aux=["dx"],
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.rk4,
        dtmax=coarse_dt,
        dt=coarse_dt,
        t_span=(0.0, 4 * pi),
        event_var="x",
        feature_var="x",
        observer_min_x_amp=0.5,
        observer_x_up_thresh=(2 + sqrt(2)) / 4,
        observer_x_down_thresh=coarse_x_down_thresh,
        observer_dx_down_thresh=0.001,
        observer_dx_up_thresh=0.001,
        observer_max_event_count=100,
        observer_max_event_timestamps=3,
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    event_count = int(output.get_var_count("event"))
    up_times = np.asarray(output.get_timestamps("up"), dtype=np.float64)
    down_times = np.asarray(output.get_timestamps("down"), dtype=np.float64)

    assert event_count == 2
    np.testing.assert_allclose(up_times, np.array([pi / 4, 9 * pi / 4]), atol=5e-2, rtol=0.0)
    x_down = -1.0 + 2.0 * coarse_x_down_thresh
    np.testing.assert_allclose(
        down_times,
        np.array([pi - np.arcsin(x_down), 3 * pi - np.arcsin(x_down)]),
        atol=6e-2,
        rtol=0.0,
    )


@pytest.mark.parametrize(
    ("event_direction", "expected_times"),
    [
        (clode.EventDirection.rising, np.array([pi / 6], dtype=np.float64)),
        (clode.EventDirection.falling, np.array([5 * pi / 6], dtype=np.float64)),
        (
            clode.EventDirection.either,
            np.array([pi / 6, 5 * pi / 6], dtype=np.float64),
        ),
    ],
)
def test_threshold_1_absolute_crossings_respect_direction(
    event_direction: clode.EventDirection,
    expected_times: np.ndarray,
) -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        observer=clode.Observer.threshold_crossing,
        stepper=clode.Stepper.rk4,
        dtmax=0.05,
        dt=0.05,
        t_span=(0.0, 2 * pi),
        event_var="x",
        feature_var="x",
        observer_event_direction=event_direction,
        observer_x_up_thresh=0.5,
        observer_max_event_count=8,
        observer_max_event_timestamps=3,
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    event_count = int(output.get_var_count("event"))
    threshold_times = np.atleast_1d(
        np.asarray(output.get_timestamps("threshold"), dtype=np.float64)
    )

    assert event_count == len(expected_times)
    np.testing.assert_allclose(threshold_times, expected_times, atol=3e-2, rtol=0.0)


def test_threshold_1_absolute_crossings_respect_min_amp_gate() -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        observer=clode.Observer.threshold_crossing,
        stepper=clode.Stepper.rk4,
        dtmax=0.05,
        dt=0.05,
        t_span=(0.0, 2 * pi),
        observer_configuration=clode.ThresholdCrossingConfig(
            event_var="x",
            threshold=0.5,
            direction=clode.EventDirection.either,
            min_amp=3.0,
            max_event_count=8,
            max_event_timestamps=3,
        ),
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    assert int(output.get_var_count("event")) == 0


@pytest.mark.parametrize(
    ("event_direction", "expected_times"),
    [
        (clode.EventDirection.rising, np.array([pi / 6], dtype=np.float64)),
        (clode.EventDirection.falling, np.array([5 * pi / 6], dtype=np.float64)),
        (
            clode.EventDirection.either,
            np.array([pi / 6, 5 * pi / 6], dtype=np.float64),
        ),
    ],
)
def test_normalized_threshold_crossing_uses_warmup_relative_thresholds(
    event_direction: clode.EventDirection,
    expected_times: np.ndarray,
) -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        observer=clode.Observer.normalized_threshold_crossing,
        observer_configuration=clode.ThresholdCrossingConfig(
            event_var="x",
            threshold=0.75,
            direction=event_direction,
            min_amp=0.5,
            max_event_count=8,
            max_event_timestamps=3,
        ),
        stepper=clode.Stepper.rk4,
        dtmax=0.05,
        dt=0.05,
        t_span=(0.0, 2 * pi),
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    event_count = int(output.get_var_count("event"))
    threshold_times = np.atleast_1d(
        np.asarray(output.get_timestamps("threshold"), dtype=np.float64)
    )

    assert event_count == len(expected_times)
    np.testing.assert_allclose(threshold_times, expected_times, atol=3e-2, rtol=0.0)


def test_normalized_threshold_crossing_respects_min_amp_gate() -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        observer=clode.Observer.normalized_threshold_crossing,
        observer_configuration=clode.ThresholdCrossingConfig(
            event_var="x",
            threshold=0.75,
            direction=clode.EventDirection.either,
            min_amp=3.0,
            max_event_count=8,
            max_event_timestamps=3,
        ),
        stepper=clode.Stepper.rk4,
        dtmax=0.05,
        dt=0.05,
        t_span=(0.0, 2 * pi),
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    assert int(output.get_var_count("event")) == 0


def test_absolute_schmitt_trigger_uses_absolute_thresholds() -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        aux=["dx"],
        observer=clode.Observer.schmitt_trigger,
        observer_configuration=clode.SchmittTriggerConfig(
            event_var="x",
            feature_var="x",
            x_up_threshold=0.5,
            x_down_threshold=-0.5,
            dx_up_threshold=0.0,
            dx_down_threshold=0.0,
            min_amp=0.5,
            max_event_count=8,
            max_event_timestamps=3,
        ),
        stepper=clode.Stepper.rk4,
        dtmax=0.05,
        dt=0.05,
        t_span=(0.0, 4 * pi),
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    up_times = np.asarray(output.get_timestamps("up"), dtype=np.float64)
    down_times = np.asarray(output.get_timestamps("down"), dtype=np.float64)
    assert int(output.get_var_count("event")) == 2
    np.testing.assert_allclose(
        up_times,
        np.array([pi / 6, 13 * pi / 6], dtype=np.float64),
        atol=3e-2,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        down_times,
        np.array([7 * pi / 6, 19 * pi / 6], dtype=np.float64),
        atol=3e-2,
        rtol=0.0,
    )


def test_local_max_coarse_timestamps_use_three_sample_refinement() -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        aux=["dx"],
        observer=clode.Observer.local_max,
        stepper=clode.Stepper.rk4,
        dtmax=0.2,
        dt=0.2,
        t_span=(0.0, 4 * pi),
        event_var="x",
        feature_var="x",
        observer_max_event_count=4,
        observer_max_event_timestamps=4,
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    max_times = np.asarray(output.get_timestamps("localmax"), dtype=np.float64)
    min_times = np.asarray(output.get_timestamps("localmin"), dtype=np.float64)

    np.testing.assert_allclose(max_times, np.array([pi / 2, 5 * pi / 2]), atol=1e-2, rtol=0.0)
    np.testing.assert_allclose(min_times, np.array([3 * pi / 2, 7 * pi / 2]), atol=1e-2, rtol=0.0)


@pytest.mark.parametrize(
    ("polarity", "expected_times", "expected_values"),
    [
        (
            clode.ExtremumPolarity.maximum,
            np.array([pi / 2, 5 * pi / 2], dtype=np.float64),
            np.array([1.0, 1.0], dtype=np.float64),
        ),
        (
            clode.ExtremumPolarity.minimum,
            np.array([3 * pi / 2, 7 * pi / 2], dtype=np.float64),
            np.array([-1.0, -1.0], dtype=np.float64),
        ),
    ],
)
def test_local_extremum_selected_polarity_uses_three_sample_refinement(
    polarity: clode.ExtremumPolarity,
    expected_times: np.ndarray,
    expected_values: np.ndarray,
) -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        aux=["dx"],
        observer=clode.Observer.local_extremum,
        observer_configuration=clode.LocalExtremumConfig(
            variable="x",
            polarity=polarity,
            max_event_count=4,
            max_event_timestamps=4,
        ),
        stepper=clode.Stepper.rk4,
        dtmax=0.2,
        dt=0.2,
        t_span=(0.0, 4 * pi),
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    event_times = np.asarray(output.get_timestamps("local_extremum"), dtype=np.float64)
    event_values = np.asarray(
        output.get_event_data("local_extremum", type="value"),
        dtype=np.float64,
    )

    assert int(output.get_var_count("event")) == 2
    np.testing.assert_allclose(event_times, expected_times, atol=1e-2, rtol=0.0)
    np.testing.assert_allclose(event_values, expected_values, atol=2e-2, rtol=0.0)
