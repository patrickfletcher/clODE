import numpy as np
import matplotlib.pyplot as plt
import clode
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
    "gk": 2.5,
    "gsk": 3,
    "gleak": 0.1,
    "cm": 4.0,
    "e_leak": -50,
    "tau_n": 15,
    "k_c": 0.1,
}

tend=8000
features_integrator = clode.FeatureSimulator(
    rhs_equation=lactotroph,
    supplementary_equations=[x_inf, s_inf],
    variables=variables,
    parameters=parameters,
    aux=['ica'],
    observer=clode.Observer.neighbourhood_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, tend),
    observer_min_x_amp=0.0,
    observer_x_down_thresh=0.02,
    observer_neighbourhood_radius=0.5,
    observer_dx_down_thresh=-0.1,
    observer_dx_up_thresh=0.1,
    observer_max_event_count=100000,
    event_var="v",
    feature_var="v",
)

features_integrator.transient()
output = features_integrator.features()

events = output.get_timestamp("event")
active = output.get_timestamp("active")
inactive = output.get_timestamp("inactive")
afterevent = output.get_timestamp("afterevent")
period_count = int(output.get_var_count("period"))
step_count = int(output.get_var_count("step"))

# Plot trajectory
trajectory_integrator = clode.TrajectorySimulator(
    rhs_equation=lactotroph,
    supplementary_equations=[x_inf, s_inf],
    variables=variables,
    parameters=parameters,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, tend),
    dt=0.01,
)

active_threshold = output.get_eventvar_threshold('x')
dxup_threshold = output.get_featurevar_threshold('dx up')
dxdown_threshold = output.get_featurevar_threshold('dx down')
inactive_v_mean = output.get_var_mean('inactive v')

trajectory_integrator.transient()
trajectory = trajectory_integrator.trajectory()

plt.plot(trajectory[0].t, trajectory[0].x[:, 0])
# Plot events
for event in events:
    plt.axvline(x=event, color="red", linestyle="--")

for event in active:
    plt.axvline(x=event, color="green", linestyle="--")

for event in inactive:
    plt.axvline(x=event, color="blue", linestyle="--")

# for event in afterevent:
#     plt.axvline(x=event, color="orange", linestyle="--")

plt.xlabel("Time")
plt.ylabel("Membrane potential")
plt.title("Lactotroph model")
plt.ylim(-80, 20)
plt.show()
pass

# Plot dv/dt by calculating the derivative of the membrane potential from trajectory.x
dvdt = np.gradient(trajectory[0].x[:, 0], trajectory[0].t)
plt.plot(trajectory[0].t, dvdt)
plt.xlabel("Time")
plt.ylabel("dv/dt")
plt.title("Lactotroph model")
for event in events:
    plt.axvline(x=event, color="red", linestyle="--")

for event in active:
    plt.axvline(x=event, color="green", linestyle="--")

for event in inactive:
    plt.axvline(x=event, color="blue", linestyle="--")
#
# for event in afterevent:
#     plt.axvline(x=event, color="orange", linestyle="--")

plt.show()
pass
