# Trajectory simulation

...Under Construction...

CLODE can simulate ODE trajectories using the CLODETrajectory class.

## Example

The following example simulates the Van der Pol oscillator using the dormand_prince45 integrator.

### OpenCL RHS function

```c
void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype m = par[0];
    realtype w = par[1];
    realtype k = par[2];
    realtype H = par[3];

    realtype y1 = var[0];
    realtype y2 = var[1];

    realtype dy1 = y2;
    realtype dy2 = (w - k * y2) / m;

    derivatives[0] = dy1;
    derivatives[1] = dy2;
    aux[0] = y1 - H;
}
```

### Python

```py run
import clode
import matplotlib.pyplot as plt
import numpy as np
from typing import List

def fitzhugh_nagumo(
    time: float,
    variables: List[float],
    parameters: List[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
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
parameters = {"a": a, "b": 0.8, "current": np.arange(0.0, 0.6, 0.1), "epsilon": 1.0 / 12.5}

simulator = clode.TrajectorySimulator(
    rhs_equation=fitzhugh_nagumo,
    variables=variables,
    parameters=parameters,
    t_span=(0, 200),
    dt=0.02,
)

trajectories = simulator.trajectory()

plt.figure(figsize=(8, 6))
for index in range(len(trajectories)):
    label = f"I={parameters['current'][index]:.1f}"
    plt.plot(trajectories[index].x[:, 0], trajectories[index].x[:, 1], label=label)
plt.xlabel("V")
plt.ylabel("w")
plt.legend()
plt.title("FitzHugh-Nagumo phase plane")
plt.show()

# Plot the time series
plt.figure(figsize=(8, 6))
for index in range(0, len(trajectories), 2):
    label = f"I={parameters['current'][index]}"
    plt.plot(trajectories[index].t, trajectories[index].x[:, 0], label=label)
plt.xlabel("t")
plt.ylabel("V")
plt.legend()
plt.title("FitzHugh-Nagumo time series")
plt.show()
```
