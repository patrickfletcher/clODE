from math import pi
from typing import List
import pytest

import clode
clode.set_log_level(clode.LogLevel.info)

# test each observer to make sure they:
# - handle zero aux
# - return expected results

# Notes:
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
# ** it takes a full period to visit the full extent of state-space (required for computing the positiong in unit-normalized state space)
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

# TODO: per-observer parameter struct?
# TODO: support zero parameters? (for sweeps of initial conditions only)

# harmonic oscillator
def get_rhs(t: float,
            vars: List[float],
            p: List[float],
            dy: List[float],
            aux: List[float],   
            w: List[float]) -> None:
    k: float = p[0]
    x: float = vars[0]
    y: float = vars[1]
    dy[0] = y
    dy[1] = -k*x

# start at the left of the unit circle
variables = {"x": -1.0, "y":0.0} 

parameters = {"k": 1.0} 
expected_period = 2*pi

# 10 full periods
expected_events = 10
t_span = (0.0, expected_events*expected_period)
dt = 0.01
expected_steps = int(t_span[1]/dt)+1

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

    integrator.features()

    print(f"final state: {integrator.get_final_state()[0]}")

    features = integrator.get_observer_results()
    feature_names = features.get_feature_names()

    #TODO: switch to asserts?
    if "mean y" in feature_names:
        mean_y = features.get_var_mean("y")
        print(f"expected mean y: {0.0},\t\t actual:{mean_y:0.4}")
    if "mean IMI" in feature_names:
        inter_maximum_interval = features.get_var_mean("IMI")
        print(f"expected IMI: {expected_period:0.4},\t\t actual:{inter_maximum_interval:0.4}")
    if "mean period" in feature_names:
        period = features.get_var_mean("period")
        print(f"expected period: {expected_period:0.4},\t\t actual:{period:0.4}")
    if "event count" in feature_names:
        event_count = features.get_var_count("event")
        print(f"expected event count: {expected_events},\t actual:{int(event_count)}")
    if "period count" in feature_names:
        period_count = features.get_var_count("period")
        print(f"expected period count: {expected_events-1},\t actual:{int(period_count)}")
    if "step count" in feature_names:
        step_count = features.get_var_count("step")
        print(f"expected step count: {expected_steps},\t actual:{int(step_count)}")
    print("\n")

observers = [observer for observer in clode.Observer]

@pytest.mark.skip(reason="for now just to validate/debug observers")
def test_all_observers():
    for observer in observers:
        test_observer(observer)

if __name__ == "__main__":
    test_all_observers()