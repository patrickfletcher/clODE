import matplotlib.pyplot as plt
import numpy as np
from typing import List

import clode
from clode import exp
# clode.set_log_level(clode.LogLevel.debug)

def get_rhs(t: float,
            x: List[float],
            p: List[float],
            dx: List[float],
            aux: List[float],   
            w: List[float]) -> None:    
    v: float = x[0]
    n: float = x[1]
    c: float = x[2]

    gca: float = p[0]
    gkca: float = p[1]
    kpmca: float = p[2]
    gk: float = 3500.0

    vca: float = 25.0
    vk: float = -75.0
    cm: float = 5300.0
    alpha: float = 4.5e-6
    fcyt: float = 0.01
    kd: float = 0.4
    vm: float = -20.0
    sm: float = 12.0
    vn: float = -16.0
    sn: float = 5.0
    taun: float = 20.0
    
    minf: float = 1.0/(1.0 + exp((vm - v)/sm))
    ninf: float = 1.0/(1.0 + exp((vn - v)/sn))
    omega: float = c**2/(c**2 + kd**2)

    ica: float = gca*minf*(v - vca)
    ik: float = gk*n*(v - vk)
    ikca: float = gkca*omega*(v - vk)

    dx[0] = -(ica + ik + ikca)/cm
    dx[1] = (ninf - n)/taun
    dx[2] = fcyt*(-alpha*ica - kpmca*c)

variables = {"v": -50.0, "n": 0.01, "c": 0.12}
parameters = {"gca": 1200.0, "gkca": 750.0, "kpmca": 0.1}

# set up the solver
integrator = clode.FeatureSimulator(
    rhs_equation=get_rhs,
    variables=variables,
    parameters=parameters,
    single_precision=True,
    stepper=clode.Stepper.dormand_prince,
    dt=0.001,
    dtmax = 0.1,
    abstol=1e-6,
    reltol=1e-5,
    event_var="v",
    feature_var="v",
    observer=clode.Observer.threshold_2,
    observer_x_up_thresh=0.5,
    observer_x_down_thresh=0.05,
)

# set up the ensemble of systems
nx = ny = 128
gca = np.linspace(550.0, 1050.0, nx)
kpmca = np.linspace(0.095, 0.155, ny)
gca_grid, kpmca_grid = np.meshgrid(gca, kpmca)

# the 2D grid shape is noted internally, and features will be returned in this shape
integrator.set_ensemble(parameters= {"gca" : gca_grid, "kpmca" : kpmca_grid})

integrator.set_tspan((0.0, 50000.0))
integrator.transient()
integrator.set_tspan((0.0, 10000.0))
integrator.features()

features = integrator.get_observer_results()

# the feature output knows to return the feature with shape matching the 2D grid input
max_peaks = features.get_var_max("peaks")

plt.pcolormesh(gca_grid, kpmca_grid, max_peaks, shading='nearest', vmax=12)
plt.title("peaks")
plt.colorbar()
plt.xlabel("gca")
plt.ylabel("kpmca")
plt.axis("tight")

# highlight a few example points - we'll get their trajectories next
points = np.array([[950, 0.145], [700, 0.105], [750, 0.125], [800, 0.142]])
plt.plot(points[:, 0], points[:, 1], 'o', color='black')
for i, txt in enumerate(range(4)):
    plt.annotate(txt, (points[i, 0] - 10, points[i, 1] - 0.003))

plt.show()


# Now get the trajectories
steps_taken = features.get_var_count("step")
max_store = int(np.max(steps_taken))

integrator_traj = clode.TrajectorySimulator(
    rhs_equation=get_rhs,
    variables = variables,
    parameters = parameters,
    single_precision = True,
    stepper = clode.Stepper.dormand_prince,
    dt = 0.001,
    dtmax = 0.1,
    abstol = 1e-6,
    reltol = 1e-5,
    max_store = max_store,
)

integrator_traj.set_ensemble(parameters = {"gca":points[:, 0], "kpmca": points[:, 1]})

integrator_traj.transient(t_span=(0.0, 50000.0))
trajectories = integrator_traj.trajectory(t_span=(0.0, 10000.0))

fig, ax = plt.subplots(4, 1, sharex=True, sharey=True)
for i, trajectory in enumerate(trajectories):
    ax[i].plot(trajectory.t / 1000.0, trajectory.x["v"])

ax[1].set_ylabel("v")
ax[-1].set_xlabel("time (s)")
plt.show()