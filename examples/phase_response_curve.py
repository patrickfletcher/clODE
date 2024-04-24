import numpy as np
import matplotlib.pyplot as plt
import clode
from clode import exp
from typing import List

# clode.set_log_level(clode.LogLevel.debug)

# 1. first pass to get the period of each oscillator
# 2. set perturbation times as linspace(0, T0), i.e. phase in [0,1]
# 3. second pass with perturbations, record k event times
# 4. post-processing to get the perturbed periods and delta phase (T0-T1)/T0 

def fitzhugh_nagumo(
    time: float,
    variables: List[float],
    parameters: List[float],
    derivatives: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    V: float = variables[0]
    w: float = variables[1]

    a: float = parameters[0]
    b: float = parameters[1]
    current: float = parameters[2]
    epsilon: float = parameters[3]

    dV: float = V - V ** 3 / 3 - w + current
    dw: float = epsilon * (V + a - b * w)

    derivatives[0] = dV
    derivatives[1] = dw

variables = {"V": 1.0, "w": 0.0}
parameters = {"a": 0.7, "b": 0.8, "current": 0.4, "epsilon": 1.0 / 12.5}

auxvars = []

tend=1000
features_integrator = clode.FeatureSimulator(
    rhs_equation=fitzhugh_nagumo,
    variables=variables,
    parameters=parameters,
    aux=auxvars,
    observer=clode.Observer.local_max,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, tend),
    observer_max_event_count=10,
    observer_max_event_timestamps = 10,
    observer_min_x_amp=0,
    observer_min_imi=0,
    feature_var="V",
    dtmax=1.0,
    dt=0.1,
    abstol=1.0e-6,
    reltol=1.0e-6
)

features_integrator.transient()
features_integrator.set_tspan((0.0, 300.))
output = features_integrator.features()

print(output)

up_times = output.get_timestamps("localmax")
down_times = output.get_timestamps("localmin")

print(up_times)
print(down_times)

# Get the trajectory
trajectory_integrator = clode.TrajectorySimulator(
    rhs_equation=fitzhugh_nagumo,
    variables=variables,
    parameters=parameters,
    aux=auxvars,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, tend),
    dtmax=1.0,
    dt=0.1,
    abstol=1.0e-6,
    reltol=1.0e-6
)

trajectory_integrator.transient()
trajectory_integrator.set_tspan((0.0, 300.))
trajectory = trajectory_integrator.trajectory()

var = "V"
t = trajectory.t
v = trajectory.x[var]
dvdt = trajectory.dx[var]

plt.plot(t, v)
# Plot events
for event in up_times:
    plt.axvline(x=event, color="red", linestyle="--")

for event in down_times:
    plt.axvline(x=event, color="blue", linestyle="--")

plt.xlabel("t")
plt.ylabel(var)
plt.show()
# pass