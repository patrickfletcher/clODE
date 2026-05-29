from typing import List

import clode
import matplotlib.pyplot as plt
import numpy as np


def van_der_pol(
    t: float,
    variables: List[float],
    parameters: List[float],
    derivatives: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    x: float = variables[0]
    y: float = variables[1]
    mu: float = parameters[0]

    derivatives[0] = y
    derivatives[1] = mu * (1.0 - x * x) * y - x


mu_values = np.array([0.01, 0.5, 2.0, 4.0])

integrator = clode.FeatureSimulator(
    rhs_equation=van_der_pol,
    variables={"x": 1.0, "y": 1.0},
    parameters={"mu": 0.1},
    observer=clode.Observer.normalized_schmitt_trigger,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
)

integrator.set_ensemble(parameters={"mu": mu_values})
integrator.transient()
observer_output = integrator.features()

period = observer_output.get_var_mean("period")
print(period)

plt.plot(mu_values, period)
plt.xlabel("mu")
plt.ylabel("period")
plt.title("Van der Pol oscillator")
plt.show()
