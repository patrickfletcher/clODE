import matplotlib.pyplot as plt
import numpy as np

import clode

src_file = "chay_keizer.cl"

variables = {"v": -50.0, "n": 0.01, "c": 0.12}
parameters = {"gca": 1200.0, "gkca": 750.0, "kpmca": 0.1}

# set up the ensemble of systems
nx = 32
ny = 32
nPts = nx * ny
gca = np.linspace(550.0, 1050.0, nx)
kpmca = np.linspace(0.095, 0.155, ny)
px, py = np.meshgrid(gca, kpmca)

ensemble_parameters = {"gca" : px.flatten(), "kpmca" : py.flatten()} #gkca will have default value
ensemble_parameters_names = list(ensemble_parameters.keys())

integrator = clode.FeatureSimulator(
    src_file=src_file,
    variables=variables,
    parameters=parameters,
    single_precision=True,
    stepper=clode.Stepper.dormand_prince,
    dt=0.001,
    dtmax=0.1,
    abstol=1e-6,
    reltol=1e-5,
    event_var="v",
    feature_var="v",
    observer=clode.Observer.threshold_2,
    observer_x_up_thresh=0.5,
    observer_x_down_thresh=0.05,
    observer_min_x_amp=5.0,
    observer_min_imi=0.0,
    observer_max_event_count=1000,
)

# X0 = np.tile(initial_state, (nPts, 1))
# X0 = X0 + (np.random.random((nPts, 1)) - 0.5) * [20.0, 0.0, 0.1]

integrator.set_ensemble(parameters=ensemble_parameters)

integrator.set_tspan((0.0, 20000.0))
integrator.transient()

integrator.set_tspan((0.0, 10000.0))
integrator.features()

features = integrator.get_observer_results()

feature = features.get_var_mean("v")
# feature = features.get_var_max("peaks")
feature = np.reshape(feature, (nx, ny))

print(feature)

plt.pcolormesh(px, py, feature, shading='nearest', vmax=12)
plt.title("peaks")
plt.colorbar()
plt.xlabel(ensemble_parameters_names[0])
plt.ylabel(ensemble_parameters_names[1])

plt.axis("tight")

points = np.array([[950, 0.145], [700, 0.105], [750, 0.125], [800, 0.142]])
plt.plot(points[:, 0], points[:, 1], 'o', color='black')

for i, txt in enumerate(range(4)):
    plt.annotate(txt, (points[i, 0] - 10, points[i, 1] - 0.003))

plt.show()


steps_taken = features.get_var_count("step")
max_steps = int(np.max(steps_taken))

print(f"max steps taken: {max_steps}")

# assert(max_steps>0)
max_steps = 20000

# Now get the trajectories

integrator_traj = clode.TrajectorySimulator(
    src_file = src_file,
    variables = variables,
    parameters = parameters,
    single_precision = True,
    stepper = clode.Stepper.dormand_prince,
    dt = 0.001,
    dtmax = 0.1,
    abstol = 1e-6,
    reltol = 1e-5,
    max_steps = max_steps + 1,  # must be int
    max_store = max_steps + 1,
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