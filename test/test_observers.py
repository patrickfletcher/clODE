from math import pi
from typing import List

import clode
clode.set_log_level(clode.LogLevel.info)

# test each observer to make sure they:
# - handle zero aux
# - return expected results

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

# start at the top of the unit circle
# - thresh2: in upstate, upThresh=0.3, downThresh=0.2
# - nhood2: should find first min   
variables = {"x": 0.0, "y":1.0} 

parameters = {"k": 1.0} 
expected_period = 2*pi

# 10 full periods
expected_events = 10
t_span = (0.0, expected_events*expected_period)
dt = 0.01
expected_steps = int(t_span[1]/dt)+1

def test_observer(observer):
    integrator = clode.FeatureSimulator(
        rhs_equation=get_rhs,
        variables=variables,
        parameters=parameters,
        dt=dt,
        t_span=t_span,
        observer=observer
    )

    integrator.set_repeat_ensemble(num_repeats=1)

    integrator.features()
    features = integrator.get_observer_results()

    feature_names = features.get_feature_names()

    #TODO: switch to asserts
    if "mean x" in feature_names:
        mean_x = features.get_var_mean("x")
        print(f"expected mean x: {0.0},\t actual:{mean_x:0.4}")
    if "mean IMI" in feature_names:
        inter_maximum_interval = features.get_var_mean("IMI")
        print(f"expected period: {expected_period:0.4},\t actual:{inter_maximum_interval:0.4}")
    if "mean period" in feature_names:
        period = features.get_var_mean("period")
        print(f"expected period: {expected_period:0.4},\t actual:{period:0.4}")
    if "event count" in feature_names:
        event_count = features.get_var_count("event")
        print(f"expected events: {expected_events},\t actual:{int(event_count)}")
    if "step count" in feature_names:
        step_count = features.get_var_count("step")
        print(f"expected steps: {expected_steps},\t actual:{int(step_count)}")
    print("\n")

observers = [observer for observer in clode.Observer]

def test_all_observers():
    for observer in observers:
        test_observer(observer)

if __name__ == "__main__":
    test_all_observers()