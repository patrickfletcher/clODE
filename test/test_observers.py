from math import cos, pi, sqrt
from typing import List

import pytest

import clode

clode.set_log_level(clode.LogLevel.info)

# test each observer to make sure they:
# - handle zero aux
# - return expected results

# oscillatory system
# harmonic oscillator. Solution is x=cos(t), y=sin(t) --> the unit circle
# Starting at the top of the unit circle, the trajectory should proceed counter-clockwise.
# We use "y" as the feature variable
#
# local_max
# - event triggers at any local maximum in event_var
# - starting at (-1.0, 0.0) means we'll see one local max in y per period
#
# nhood1 (one-pass, neighborhood event detector) **first period can be off
# - trigger point is first local min of event_var
# - after that, event triggers if the trajectory enters a ball of radius=observer_neighbourhood_radius, measured in unit-normalized state space
# - starting at (-1.0, 0.0) should trigger on first rotation.
# ** it takes a full period to visit the full extent of state-space (required for computing the position in unit-normalized state space)
# ** this one doesn't do a warmup pass, so the unit-normalization will be based on unpredictable amount of event_var range! here, exactly half of truth...
# **--> however, it adapts: as long as it runs a full period, the threshold will settle to correct value.
# - number of periods = event count - 1
#
# nhood2 (two-pass, neighborhood event detector)
# - trigger point is the point in state space where eVarIx first drops below observer_x_down_thresh
# - after that, event triggers if the trajectory enters a ball of radius=observer_neighbourhood_radius, measured in unit-normalized state space
# - a warmup pass is done to measure the full extent of state-space, required for computing the positiong in unit-normalized state space
# - starting at (-1.0, 0.0) should trigger on first rotation.
# - number of periods = event count - 1
#
# thresh2 (two-pass, Shmitt trigger-like threshold detector)
# - thresholds calculated based on unit-normalized values of event_var (first pass finds max and min for normalization)
# - initialization: if event_var > observer_x_up_thresh, considered to be in the "up-state"
# - event triggers at upward crossing of observer_x_up_thresh for the event_var (from "down-state" to "up-state")
# - starting at (-1.0, 0.0) should start in the up-state; triggers on first rotation
# - number of periods = event count - 1

# harmonic oscillator. 
# Solution: x=cos(t), y=sin(t)
def get_rhs(
    t: float,
    vars: List[float],
    p: List[float],
    dy: List[float],
    aux: List[float],
    w: List[float],
) -> None:
    k: float = p[0]
    x: float = vars[0]
    y: float = vars[1]
    dy[0] = y
    dy[1] = -k * x


# start at the left of the unit circle
variables = {"x": -1.0, "y": 0.0}

parameters = {"k": 1.0}
expected_period = 2 * pi

# 10 full periods
expected_events = 10
t_span = (0.0, expected_events * expected_period)
dt = 0.01
expected_steps = int(t_span[1] / dt) + 1


@pytest.mark.skip(reason="for now just to validate/debug observers")
def test_observer(observer):
    integrator = clode.FeatureSimulator(
        rhs_equation=get_rhs,
        variables=variables,
        parameters=parameters,
        dt=dt,
        t_span=t_span,
        observer=observer,
        event_var="y",
        feature_var="y",
        observer_min_x_amp=0.0,
        observer_x_up_thresh=0.3,
        observer_x_down_thresh=0.2,
        observer_neighbourhood_radius=0.05,
    )

    integrator.set_repeat_ensemble(num_repeats=1)

    features = integrator.features()
    feature_names = features.get_feature_names()

    print(f"--- {observer} ---")
    print(f"final state: {integrator.get_final_state()[0]}")

    # TODO: map of feature_names --> expected results
    # TODO: switch to asserts?
    if "mean x" in feature_names:
        mean_x = features.get_var_mean("x")
        print(f"expected mean x: {0.0},\t\t actual:{mean_x:0.6}")
    if "max x" in feature_names:
        max_x = features.get_var_max("x")
        print(f"expected max x: {1.0},\t\t actual:{max_x:0.6}")
    if "min x" in feature_names:
        min_x = features.get_var_min("x")
        print(f"expected min x: {-1.0},\t\t actual:{min_x:0.6}")
    if "mean y" in feature_names:
        mean_y = features.get_var_mean("y")
        print(f"expected mean y: {0.0},\t\t actual:{mean_y:0.6}")
    if "mean amplitude" in feature_names:
        amp = features.get_var_mean("amplitude")
        print(f"expected amplitude: {2.0},\t actual:{amp:0.6}")
    if "mean IMI" in feature_names:
        imi = features.get_var_mean("IMI")
        print(f"expected IMI: {expected_period:0.6},\t\t actual:{imi:0.6}")
    if "mean period" in feature_names:
        period = features.get_var_mean("period")
        print(f"expected period: {expected_period:0.6},\t actual:{period:0.6}")
    if "event count" in feature_names:
        event_count = features.get_var_count("event")
        print(f"expected event count: {expected_events},\t actual:{int(event_count)}")
    if "step count" in feature_names:
        step_count = features.get_var_count("step")
        print(f"expected step count: {expected_steps},\t actual:{int(step_count)}")
    print("\n")


@pytest.mark.skip(reason="for now just to validate/debug observers")
def test_all_observers():
    observers = [observer for observer in clode.Observer]
    for observer in observers:
        test_observer(observer)


def test_threshold_2_observes_sine_events():
    # Follow clode rhs convention
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

    # Define the parameters
    parameters = {
        "dilation": 1,
    }

    # Define the initial conditions
    variables = {
        "x": 0,
    }

    # Activate at t=pi/4, deactivate at t=3pi/2
    features_integrator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables=variables,
        parameters=parameters,
        observer=clode.Observer.threshold_2,
        aux=["dx"],
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 4 * pi),
        observer_min_x_amp=0.5,
        observer_x_up_thresh=(2 + sqrt(2)) / 4,
        observer_x_down_thresh=0.001,
        observer_dx_down_thresh=0.001,
        observer_dx_up_thresh=0.001,
        event_var="x",
        feature_var="x",
        dtmax=0.001,
        dt=0.001,
    )

    # Run the simulation
    output = features_integrator.features()

    # Get the number of periods (should be 1, as one event has not finished)
    event_count = int(output.get_var_count("period"))
    assert event_count == 1

    # Get the timestamps of the events
    up_times = output.get_timestamps("up")
    assert up_times[0] == pytest.approx(pi / 4, rel=1e-2)
    assert up_times[1] == pytest.approx(9 * pi / 4, rel=1e-2)
    down_times = output.get_timestamps("down")
    assert down_times[0] == pytest.approx(3 * pi / 2, rel=1e-2)
    assert down_times[1] == pytest.approx(7 * pi / 2, rel=1e-2)


if __name__ == "__main__":
    test_all_observers()
