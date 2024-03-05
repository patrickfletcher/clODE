# Stochastic simulations

Stochastic simulations are supported via the Euler-Maruyama method. Random normal variables $W_t\sim N(0, dt)$ can be included in the vector field function, and each simulation instance will have an independent stream of random variables.

As an example, consider the [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein-Uhlenbeck_process) with $\theta = 1$:

$$ dx_t = (\mu - x_t)dt + \sigma W_t $$

To implement this system, we heuristically write it in Langevin form:

$$ \frac{dx}{dt} = \mu - x + \sigma\eta(t)$$

where $\eta(t)$ represents white noise. When $\sigma=0$, the system is a first-order ODE with steady state x=$\mu$. Increasing $\sigma$ produces solutions with a steady state distribution $x\sim N(\mu, \frac{\sigma^2}{2})$.

```py run
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

t_span = (0.0, 1000.0)

integrator = clode.Simulator(
    rhs_equation=wiener_process,
    variables=variables,
    parameters=parameters,
    num_noise=1,
    stepper=clode.Stepper.stochastic_euler,
    single_precision=True,
    t_span=t_span,
    dt=0.001,
)

# set up the ensemble of Wiener processes with identical parameters and initial conditions
# Any parameters or initial conditions that are not specified will be set to the default values
nPts = 8192
integrator.set_repeat_ensemble(nPts)

integrator.transient()

XF = integrator.get_final_state()

plt.hist(XF, 30)
plt.xlabel("x")
plt.show()

print(f"mean xf: {np.mean(XF) :0.5}")
print(f"simulation variance: {np.var(XF) :0.5}")
print(f"expected variance: {parameters['sigma'] ** 2 / 2 :0.5}")
```

Output:

```bash
mean xf: 1.0041
simulation variance: 0.12629
expected variance: 0.125
```

[//]: # (![Result]&#40;Ornstein-Uhlenbeck_process.png&#41;)
