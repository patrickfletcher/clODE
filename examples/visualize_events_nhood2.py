import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
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
    "v": -60.0,
    "n": 0.1,
    "c": 0.1,
}

parameters = {
    "gca": 2.5,
    "gk": 2.56,
    "gsk": 3.0,
    "gleak": 0.1,
    "cm": 4.0,
    "e_leak": -50.0,
    "tau_n": 15.0,
    "k_c": 0.1,
}

auxvars = ['ica']

tend=2000
nhood_radius=0.05
features_integrator = clode.FeatureSimulator(
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
    observer=clode.Observer.neighbourhood_2,
    event_var="c",
    feature_var="v",
    observer_x_down_thresh=0.05,
    observer_neighbourhood_radius=nhood_radius,
    observer_max_event_count=10,
    observer_max_event_timestamps = 10,
)

features_integrator.transient()
output = features_integrator.features()

event_times = output.get_event_data("nhood","time")

print(event_times)

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
    dtmax=2.0,
    abstol=1.0e-6,
    reltol=1.0e-4,
)

trajectory_integrator.transient()
trajectory = trajectory_integrator.trajectory()

vars = list(variables.keys())

t = trajectory.t
v = trajectory.x["v"]
n = trajectory.x["n"]
c = trajectory.x["c"]

# trajectory 
plt.plot(t, v)
# Plot events
for event in event_times:
    plt.axvline(x=event, color="red", linestyle="--")

plt.xlabel("t")
plt.ylabel("v")
plt.show()

# Parametric equations for unit sphere, scale each axis and translate to nhood center
center_v = output.F["nhood center v"]
center_n = output.F["nhood center n"]
center_c = output.F["nhood center c"]
range_v = output.F["range v"]
range_n = output.F["range n"]
range_c = output.F["range c"]
phi = np.linspace(0, np.pi, 100)
theta = np.linspace(0, 2*np.pi, 100)
phi, theta = np.meshgrid(phi, theta)
x = center_c + range_c * nhood_radius * np.sin(phi) * np.cos(theta)
y = center_n + range_n * nhood_radius * np.sin(phi) * np.sin(theta)
z = center_v + range_v * nhood_radius * np.cos(phi)

# Plot trajectory with event markers in state space
event0_idx = np.argmax(t>=event_times[0])
event1_idx = np.argmax(t>=event_times[1])

ax = plt.subplot(1, 1, 1, projection='3d')
ax.plot(c[event0_idx-1:event1_idx], n[event0_idx-1:event1_idx], v[event0_idx-1:event1_idx],'.-')
ax.plot(c[event0_idx], n[event0_idx], v[event0_idx],'r>')
ax.plot(center_c, center_n, center_v, 'ko')
ax.plot_surface(x, y, z, color='linen', alpha=0.5)

ax.set_xlabel("c")
ax.set_ylabel("n")
ax.set_zlabel("v")
plt.show()
# pass