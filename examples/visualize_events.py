import numpy as np
import matplotlib.pyplot as plt
import clode
from clode import exp
from typing import List

def x_inf(v: float, vx: float, sx: float) -> float:
    return 1.0 / (1.0 + exp((vx - v) / sx))

def s_inf(c: float, k_s: float) -> float:
    c2: float = c * c
    return c2 / (c2 + k_s * k_s)

def lactotroph(
    t: float,
    x_: List[float],
    p_: List[float],
    dx_: List[float],
    aux_: List[float],
    w_: List[float],
) -> None:
    v: float = x_[0]
    n: float = x_[1]
    c: float = x_[2]

    gca: float = p_[0]
    gk: float = p_[1]
    gsk: float = p_[2]
    gleak: float = p_[3]
    cm: float = p_[4]
    e_leak: float = p_[5]
    tau_n: float = p_[6]
    k_c: float = p_[7]

    e_ca: float = 60
    e_k: float = -75

    vm: float = -20
    vn: float = -5
    sm: float = 12
    sn: float = 10

    f_c: float = 0.01
    alpha: float = 0.0015
    k_s: float = 0.4

    ica: float = gca * x_inf(v, vm, sm) * (v - e_ca)
    ik: float = gk * n * (v - e_k)
    isk: float = gsk * s_inf(c, k_s) * (v - e_k)
    ileak: float = gleak * (v - e_leak)
    current: float = ica + ik + isk + ileak

    dv: float = -current / cm
    dn: float = (x_inf(v, vn, sn) - n) / tau_n
    dc: float = -f_c * (alpha * ica + k_c * c)

    dx_[0] = dv
    dx_[1] = dn
    dx_[2] = dc
    aux_[0] = ica
    # aux_[1] = ik

clode.set_log_level(clode.LogLevel.warn)

variables = {
    "v": -60,
    "n": 0.1,
    "c": 0.1,
}

parameters = {
    "gca": 2.5,
    "gk": 2.56,
    "gsk": 3,
    "gleak": 0.1,
    "cm": 4.0,
    "e_leak": -50,
    "tau_n": 15,
    "k_c": 0.1,
}

auxvars = ['ica']

tend=4500
features_integrator = clode.FeatureSimulator(
    rhs_equation=lactotroph,
    supplementary_equations=[x_inf, s_inf],
    variables=variables,
    parameters=parameters,
    aux=auxvars,
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, tend),
    observer_min_x_amp=0.1,
    observer_x_up_thresh=0.2,
    observer_dx_up_thresh=0.,
    observer_x_down_thresh=0.1,
    observer_dx_down_thresh=0.,
    observer_max_event_count=2,
    event_var="v",
    feature_var="v",
    dtmax=1.0,
    dt=0.1,
    abstol=1.0e-6,
    reltol=1.0e-6
)

features_integrator.transient()
output = features_integrator.features()

# events = output.get_timestamps("event")
up_times = output.get_timestamps("up")
down_times = output.get_timestamps("down")
period_count = int(output.get_var_count("period"))
step_count = int(output.get_var_count("step"))
active_dip = output.get_var_mean("activeDip")

print(up_times)

# Get the trajectory
trajectory_integrator = clode.TrajectorySimulator(
    rhs_equation=lactotroph,
    supplementary_equations=[x_inf, s_inf],
    variables=variables,
    parameters=parameters,
    aux=auxvars,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, tend),
    dt=0.1,
    dtmax=1.0,
    abstol=1.0e-6,
    reltol=1.0e-6
)

trajectory_integrator.transient()
trajectory = trajectory_integrator.trajectory()

print(trajectory)

var = "v"
t = trajectory[0].t
v = trajectory[0].x(var)
dvdt = trajectory[0].dx(var)

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


# plot the dvdt vs v for the first full period
start = np.argmax(t>=up_times[0])
stop = np.argmax(t>=up_times[1])

plt.plot(v[start:stop], dvdt[start:stop])
plt.xlabel(var)
plt.ylabel(f"d{var}/dt")
plt.title("Lactotroph model")
plt.axhline(y=0.0, color="gray", linestyle="-")

# this part won't be needed if we return the state/slope at each event
first_active_idx = np.argmax(t>up_times[0]) 
first_inactive_idx = np.argmax(t>down_times[0]) 

plt.axvline(x=v[first_active_idx], color="red", linestyle="--")
plt.axhline(y=dvdt[first_active_idx], color="orange", linestyle="--")

plt.axvline(x=v[first_inactive_idx], color="blue", linestyle="--")
plt.axhline(y=dvdt[first_inactive_idx], color="green", linestyle="--")

plt.show()

# from functools import reduce
# start = 221
# end = 352
# zz = lambda: zip(trajectory[0].x[start:end, 0], trajectory[0].t[start:end])
# t0 = trajectory[0].t[start]
# mvn_avg = reduce(lambda x, y: (x[0] + (y[0] - x[0]) / (y[1] - t0), y[1]), zz())
# mvn_avg2 = reduce(lambda x, y: (x[0] + (y[0] - x[0]) * (y[1] - x[1]) / (y[1] - t0), y[1]), zz())


# pass
