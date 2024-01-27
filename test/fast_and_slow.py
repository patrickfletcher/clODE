from typing import List

import matplotlib.pyplot as plt
import numpy as np

import clode

# The FitzHugh-Nagumo model


def fitzhugh_nagumo(
    time: float,
    variables: List[float],
    parameters: List[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
) -> None:
    v: float = variables[0]
    w: float = variables[1]

    a: float = parameters[0]
    b: float = parameters[1]
    current: float = parameters[2]
    epsilon: float = parameters[3]

    dv: float = v - v**3 / 3 - w + current
    dw: float = epsilon * (v + a - b * w)

    derivatives[0] = dv
    derivatives[1] = dw


# a = [0.7, 0.8, 0.9, 1.0] * 5
#
# simulator = clode.FeatureSimulator(
#     rhs_equation=fitzhugh_nagumo,
#     variables={"v": np.arange(-2, 2, 0.2), "w": -1.0},
#     parameters={"a": a, "b": 0.8, "epsilon": 0.01, "current": 1.0},
#     observer=clode.Observer.threshold_2
# )
#
# simulator.transient()
#
# observer_output = simulator.features(True)
#
# amplitude = observer_output.get_var_max("v") - observer_output.get_var_min("v")
# period_min = observer_output.get_var_min("period")
#
# # Object classes, silent=0, excitable=1, spiking=2
# classes = np.zeros(len(amplitude))
#
# classes[amplitude > 20] = 1  # excitable
# classes[(amplitude > 20) & (period_min > 5)] = 2  # spiking
#
# # Plotting
# plt.figure(figsize=(8, 6))
# plt.scatter(a, voltages, c=classes, cmap="viridis")
#
# plt.show()
#
# print("Foo")
#
#
# a = [0.7, 0.8, 0.9, 1.0] * 5
#
# simulator = clode.TrajectorySimulator(
#     rhs_equation=fitzhugh_nagumo,
#     variables={"v": np.arange(-2, 2, 0.2), "w": -1.0},
#     parameters={"a": a, "b": 0.8, "epsilon": 0.01, "current": 1.0},
# )
#
# simulator.transient()
#
# trajectories = simulator.trajectory()
# print("Bar")
