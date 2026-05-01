import numpy as np
import matplotlib.pyplot as plt
import clode
from typing import List


# Define the Wiener_process
def wiener_process(
        t: float,
        y: List[float],
        p: List[float],
        dy: List[float],
        aux: List[float],
        wiener: List[float],
) -> None:
    x: float = y[0]
    mu: float = p[0]
    sigma: float = p[1]
    weiner_variable: float = wiener[0]
    dx: float = mu - x + sigma * weiner_variable
    dy[0] = dx


variables = {"x": 0.0}
parameters = {"mu": 1.0, "sigma": 0.5}

t_span = (0.0, 20.0)

integrator = clode.Simulator(
    rhs_equation=wiener_process,
    variables=variables,
    parameters=parameters,
    num_noise=1,
    stepper=clode.Stepper.stochastic_euler,
    single_precision=True,
    t_span=t_span,
    dt=0.01,
)

# set up the ensemble of Wiener processes with identical parameters and initial state
# Any parameters or initial state that are not specified will be set to the default values
nPts = 8192
integrator.set_repeat_ensemble(nPts)

integrator.transient()

XF = integrator.get_final_state()

print(f"expected mean: {parameters['mu'] :0.5}, simulation mean: {np.mean(XF) :0.5}")
print(f"expected variance: {parameters['sigma'] ** 2 / 2 :0.5}, simulation variance: {np.var(XF) :0.5}")

plt.hist(XF, 30)
plt.xlabel("x")
plt.show()