import numpy as np
import matplotlib.pyplot as plt
import clode
from clode import exp, heaviside
from typing import List

# clode.set_log_level(clode.LogLevel.debug)

# 1. first pass to get the period of each oscillator
# 2. set perturbation times as linspace(0, T0), i.e. phase in [0,1]
# 3. second pass with perturbations, record k event times
# 4. post-processing to get the perturbed periods and delta phase (T0-T1)/T0 

def fitzhugh_nagumo(
    t: float,
    x: List[float],
    p: List[float],
    dx: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    v: float = x[0]
    w: float = x[1]

    a: float = 0.7
    b: float = 0.8
    epsilon: float = 1.0 / 12.5

    strength: float = p[0]
    t_on: float = p[1]
    duration: float = p[2]

    t_off: float = t_on+duration
    stim: float = strength * (heaviside(t-t_on)-heaviside(t-t_off))

    dv: float = v - v ** 3 / 3 - w + 0.4 + stim
    dw: float = epsilon * (v + a - b * w)

    dx[0] = dv
    dx[1] = dw

variables = {"v": 1.0, "w": 0.0}
parameters = {"strength": 0., "t_on": 0.0, "duration": 1.0}
auxvars = []

# Set up the solver. 
tend=2000
feature_simulator = clode.FeatureSimulator(
    rhs_equation=fitzhugh_nagumo,
    variables=variables,
    parameters=parameters,
    aux=auxvars,
    stepper=clode.Stepper.rk4,
    t_span=(0.0, tend),
    dt=0.01,
    observer=clode.Observer.neighbourhood_2,
    observer_max_event_count = 4,
    observer_max_event_timestamps = 4,
    feature_var="v",
)

# Integrate first to forget the initial condition
feature_simulator.transient()

# Now we integrate using events to stop the simulation at the beginning of a period
output = feature_simulator.features()

period = output.F["mean period"][0]
print(period)

# make the grid of stimulus inputs: vary t_on, strength(?)
num_perturbation_times = 128
on_times = period + np.linspace(0, period, num_perturbation_times)
strength = 0.08

# the final state from previous simulation can automatically be used for the new ensemble initial condition
feature_simulator.set_ensemble(parameters = {"t_on": on_times, "strength": strength})

output = feature_simulator.features(t_span=(0.0, 5*period))

event_times = output.get_event_data("nhood","time")

periods_perturb = np.diff(event_times, axis=1)
stim_phase = (on_times - period)/period 
delta_phase = (period - periods_perturb)/period

ax = plt.subplot(1, 1, 1)
ax.plot(stim_phase, delta_phase[:,0])
plt.xlabel("phase ($\phi$)")
plt.ylabel("$\Delta\phi$")
plt.show()