from typing import List

import matplotlib.pyplot as plt
import numpy as np

import clode


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

    dV: float = V - V**3 / 3 - w + current
    dw: float = epsilon * (V + a - b * w)

    derivatives[0] = dV
    derivatives[1] = dw


a = 0.7
variables = {"V": 1.0, "w": 0.0}
parameters = {
    "a": a,
    "b": 0.8,
    "current": 0.0,
    "epsilon": 1.0 / 12.5,
}

simulator = clode.TrajectorySimulator(
    rhs_equation=fitzhugh_nagumo,
    variables=variables,
    parameters=parameters,
    t_span=(0, 200),
    dt=0.02,
)

ensemble_parameters = {"current":np.arange(0.0, 0.6, 0.1)}
simulator.set_ensemble(parameters=ensemble_parameters)

trajectories = simulator.trajectory()

plt.figure(figsize=(8, 6))
for index in range(len(trajectories)):
    plt.plot(
        trajectories[index].x[:, 0],
        trajectories[index].x[:, 1],
        label=f"I={ensemble_parameters['current'][index]}",
    )
plt.xlabel("V")
plt.ylabel("w")
plt.legend()
plt.show()

# Plot the time series
plt.figure(figsize=(8, 6))
for index in range(0, len(trajectories), 2):
    plt.plot(
        trajectories[index].t,
        trajectories[index].x[:, 0],
        label=f"I={ensemble_parameters['current'][index]}",
    )
plt.xlabel("t")
plt.ylabel("V")
plt.legend()
plt.show()

# import clode
# from clode import OpenCLConverter, OpenCLExp
#
#
# def helper_function(x: float, y: float) -> float:
#     return (x - y) / OpenCLExp(y)
#
#
# def get_rhs(
#     t: float,
#     variables: list[float],
#     parameters: list[float],
#     derivatives: list[float],
#     aux: list[float],
#     wiener: list[float],
# ) -> None:
#     x: float = variables[0]
#     y: float = variables[1]
#
#     derivatives[0] = helper_function(x, y)
#     derivatives[1] = -x
#
#
# converter = OpenCLConverter()
# _ = converter.convert_to_opencl(helper_function)
# opencl = converter.convert_to_opencl(get_rhs)
#
# simulator = clode.TrajectorySimulator(
#     rhs_equation=get_rhs,
#     variables={"x": 0.0, "y": 2.0},
#     parameters={"a": 1.0},
#     supplementary_equations=[helper_function],
#     t_span=(0, 10),
#     dt=0.01,
# )
#
# trajs = simulator.trajectory()
#
#
# print(opencl)
