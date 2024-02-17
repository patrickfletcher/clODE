import matplotlib.pyplot as plt
import numpy as np
from typing import List

import clode
from clode import exp

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
    observer_min_x_amp=1.0,
    observer_min_imi=0.0,
    observer_max_event_count=50,
)

# set up the ensemble of systems
nx = 512
ny = 512
nPts = nx * ny
gca = np.linspace(550.0, 1050.0, nx)
kpmca = np.linspace(0.095, 0.155, ny)
px, py = np.meshgrid(gca, kpmca)

ensemble_parameters = {"gca" : px.flatten(), "kpmca" : py.flatten()} #gkca will have default value
ensemble_parameters_names = list(ensemble_parameters.keys())

integrator.set_ensemble(parameters=ensemble_parameters)

integrator.set_tspan((0.0, 50000.0))
integrator.transient()
integrator.set_tspan((0.0, 10000.0))
integrator.features()

features = integrator.get_observer_results()

feature = features.get_var_max("peaks")
feature = np.reshape(feature, (nx, ny))

plt.pcolormesh(px, py, feature, shading='nearest', vmax=12)
plt.title("peaks")
plt.colorbar()
plt.xlabel(ensemble_parameters_names[0])
plt.ylabel(ensemble_parameters_names[1])
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

traj_parameters = {"gca":points[:, 0], "kpmca": points[:, 1]}

integrator_traj.set_ensemble(parameters = traj_parameters)

integrator_traj.set_tspan((0.0, 50000.0))
integrator_traj.transient()
integrator_traj.set_tspan((0.0, 10000.0))
integrator_traj.trajectory()

trajectories = integrator_traj.get_trajectory()

fig, ax = plt.subplots(4, 1, sharex=True, sharey=True)
for i, trajectory in enumerate(trajectories):
    ax[i].plot(trajectory.t / 1000.0, trajectory.x[:, 0])

ax[1].set_ylabel("v")
ax[-1].set_xlabel('time (s)')
plt.show()