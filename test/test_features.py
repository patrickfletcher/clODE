from math import cos, exp, pi, sqrt
from typing import List

import numpy as np

import clode


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

    # Activate at t=pi/4, deactivate at t=3pi/2
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

    assert np.isclose(up_times[0], pi / 4, atol=0.01)
    assert np.isclose(up_times[1], 9 * pi / 4, atol=0.01)

    assert np.isclose(down_times[0], 3 / 2 * pi, atol=0.01)
    assert np.isclose(down_times[1], 7 / 2 * pi, atol=0.01)

    up_times_timestamps = output.get_timestamps("up")
    down_times_timestamps = output.get_timestamps("down")

    assert len(up_times_timestamps) == len(up_times)
    assert len(down_times_timestamps) == len(down_times)

    for index in range(len(up_times)):
        assert up_times_timestamps[index] == up_times[index]

    for index in range(len(down_times)):
        assert down_times_timestamps[index] == down_times[index]
    pass
