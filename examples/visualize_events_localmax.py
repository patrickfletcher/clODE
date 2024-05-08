import numpy as np
import matplotlib.pyplot as plt
import clode
from clode import exp
from typing import List

# clode.set_log_level(clode.LogLevel.debug)

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
    "v": -60.0,
    "n": 0.1,
    "c": 0.1,
}

parameters = {
    "gca": 2.0,
    "gk": 3.0,
    "gsk": 1.5,
    "gleak": 0.,
    "cm": 4.0,
    "e_leak": -50.0,
    "tau_n": 30.0,
    "k_c": 0.1,
}

auxvars = ['ica']

tend=2000
feature_simulator = clode.FeatureSimulator(
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
    reltol=1.0e-4,
    observer=clode.Observer.local_max,
    feature_var="v",
    observer_max_event_count=10,
    observer_max_event_timestamps = 10,
)

feature_simulator.transient()
output = feature_simulator.features()

localmax_times = output.get_event_data("localmax","time")
localmax_evars = output.get_event_data("localmax","evar")
localmin_times = output.get_event_data("localmin","time")
localmin_evars = output.get_event_data("localmin","evar")

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
    reltol=1.0e-4,
)

trajectory_integrator.transient()
trajectory = trajectory_integrator.trajectory()

var = "v"
t = trajectory.t
v = trajectory.x[var]
dvdt = trajectory.dx[var]

# Plot trajectory with event markers
ax = plt.subplot(1, 1, 1)
ax.plot(t, v)
ax.plot(localmax_times, localmax_evars, 'rv')
ax.plot(localmin_times, localmin_evars, 'b^')

plt.xlabel("t")
plt.ylabel(var)
plt.show()
# pass