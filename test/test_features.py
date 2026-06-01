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


def dual_sine_curve(
    t: float,
    x_: List[float],
    p_: List[float],
    dx_: List[float],
    aux_: List[float],
    w_: List[float],
) -> None:
    dilation: float = p_[0]
    base_dx: float = cos(t * dilation)
    dx_[0] = base_dx
    dx_[1] = 2.0 * base_dx



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
        np.asarray(output.get_timestamps("event"), dtype=np.float64)
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
        np.asarray(output.get_timestamps("event"), dtype=np.float64)
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


def test_threshold_crossing_exposes_event_times_and_event_count() -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
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
        dtmax=0.05,
        dt=0.05,
        t_span=(0.0, 2 * pi),
        **device_kwargs_for_tests(),
    )

    output = feature_simulator.features()

    assert output is not None
    assert output.get_feature_names() == [
        "event time 0",
        "event time 1",
        "event time 2",
        "event count",
        "period max",
        "period min",
        "period mean",
        "n maxima max",
        "n maxima min",
        "n maxima mean",
        "amplitude max",
        "amplitude min",
        "amplitude mean",
        "max x",
        "min x",
        "mean x",
        "max dx/dt",
        "min dx/dt",
    ]
    assert int(output.get_var_count("event")) == 2


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
    assert output.get_feature_names() == [
        "up transition time 0",
        "down transition time 0",
        "up transition time 1",
        "down transition time 1",
        "up transition time 2",
        "down transition time 2",
        "event count",
        "period max",
        "period min",
        "period mean",
        "n maxima max",
        "n maxima min",
        "n maxima mean",
        "up duration max",
        "up duration min",
        "up duration mean",
        "down duration max",
        "down duration min",
        "down duration mean",
        "duty max",
        "duty min",
        "duty mean",
        "active dip max",
        "active dip min",
        "active dip mean",
        "amplitude max",
        "amplitude min",
        "amplitude mean",
        "max x",
        "min x",
        "mean x",
        "max dx/dt",
        "min dx/dt",
        "max dx",
        "min dx",
        "mean dx",
    ]
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
    np.testing.assert_allclose(output.get_var_mean("period"), 2 * pi, atol=5e-2, rtol=0.0)
    np.testing.assert_allclose(output.get_var_mean("up duration"), pi, atol=5e-2, rtol=0.0)
    np.testing.assert_allclose(output.get_var_mean("down duration"), pi, atol=5e-2, rtol=0.0)
    np.testing.assert_allclose(output.get_var_mean("duty"), 0.5, atol=2e-2, rtol=0.0)


def test_schmitt_trigger_measures_amplitude_on_feature_var() -> None:
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=dual_sine_curve,
        variables={"x": 0.0, "y": 0.0},
        parameters={"dilation": 1.0},
        observer=clode.Observer.schmitt_trigger,
        observer_configuration=clode.SchmittTriggerConfig(
            event_var="x",
            feature_var="y",
            x_up_threshold=0.5,
            x_down_threshold=-0.5,
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
    np.testing.assert_allclose(output.get_var_mean("amplitude"), 4.0, atol=1e-1, rtol=0.0)


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
    max_times = np.asarray(output.get_timestamps("local maximum"), dtype=np.float64)
    min_times = np.asarray(output.get_timestamps("local minimum"), dtype=np.float64)

    np.testing.assert_allclose(max_times, np.array([pi / 2, 5 * pi / 2]), atol=1e-2, rtol=0.0)
    np.testing.assert_allclose(min_times, np.array([3 * pi / 2, 7 * pi / 2]), atol=1e-2, rtol=0.0)


def test_local_maximum_detects_maxima_with_three_sample_refinement() -> None:
    """Verify that local_max detects local maxima with quadratic refinement."""
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables={"x": 0.0},
        parameters={"dilation": 1.0},
        aux=["dx"],
        observer=clode.Observer.local_max,
        observer_configuration=clode.LocalMaximumConfig(
            event_var="x",
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
    event_times = np.asarray(output.get_timestamps("local maximum"), dtype=np.float64)
    event_values = np.asarray(
        output.get_event_data("local maximum", type="value"),
        dtype=np.float64,
    )

    # local_max should detect maxima at pi/2 and 5*pi/2 for sine.
    expected_times = np.array([pi / 2, 5 * pi / 2], dtype=np.float64)
    expected_values = np.array([1.0, 1.0], dtype=np.float64)
    
    assert int(output.get_var_count("event")) == 2
    np.testing.assert_allclose(event_times, expected_times, atol=1e-2, rtol=0.0)
    np.testing.assert_allclose(event_values, expected_values, atol=2e-2, rtol=0.0)
