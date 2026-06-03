# Getting Started

## Basic concepts

clODE solves ensembles of ordinary differential equations (ODEs) using OpenCL.
The main user-facing classes are:

- `Simulator`: advance an ensemble and keep only the final state.
- `FeatureSimulator`: compute trajectory features on the fly through a stateful observer, without storing the full trajectory.
- `TrajectorySimulator`: store full trajectories for later inspection or plotting.

The main workflow is:

1. define an ODE as a fully typed Python function or load a `.cl` or `.xpp` source file
2. construct a simulator with default variables and parameters
3. broadcast an ensemble with `set_ensemble(...)`
4. run `transient()`, `features()`, or `trajectory()`

Repeated runs with the default `update_x0=True` continue both state and the nominal
time window. See [continuation.md](continuation.md) when you want to reset the
timebase, branch from a saved state, or inspect the difference between
`get_tspan()` and `get_final_time()`.

## Example: Van der Pol summaries without stored trajectories

[The Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator) is

$$
\dot{x} = y, \qquad \dot{y} = \mu (1 - x^2) y - x
$$

In clODE, a Python RHS must be fully typed and must use the OpenCL-converter signature shown below.

```python
from typing import List


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
```

The first example below uses the default summary observer to measure one simple readout across an ensemble of `mu` values without storing full trajectories. Richer event observers for periods, threshold crossings, and extrema are covered in [feature_extraction.md](feature_extraction.md).

```python
from typing import List

import clode
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


t_span = (0.0, 1000.0)
mu_values = np.array([0.01, 0.5, 2.0, 4.0])

integrator = clode.FeatureSimulator(
    rhs_equation=van_der_pol,
    variables={"x": 1.0, "y": 1.0},
    parameters={"mu": 0.1},
    t_span=t_span,
)

integrator.set_ensemble(parameters={"mu": mu_values})
summary_output = integrator.features()
print(summary_output.get_var_max("x"))
```

Because `FeatureSimulator` defaults to `Observer.summary`, this is the lightest-weight way to compute means, minima, and maxima on the device. When you want periods, threshold events, or local maxima, continue to [feature_extraction.md](feature_extraction.md).

## Trajectories

Use `TrajectorySimulator` when you want time-series samples rather than observer reductions.

```python
import clode
import matplotlib.pyplot as plt
import numpy as np


trajectory_simulator = clode.TrajectorySimulator(
    rhs_equation=van_der_pol,
    variables={"x": 1.0, "y": 1.0},
    parameters={"mu": 0.1},
    t_span=t_span,
)

trajectory_simulator.set_ensemble(parameters={"mu": mu_values})
trajectory_simulator.transient()
trajectories = trajectory_simulator.trajectory()

for trajectory in trajectories:
    plt.plot(trajectory.t, trajectory.x["x"])

plt.xlabel("time")
plt.ylabel("x")
plt.show()
```

Each `TrajectoryOutput` exposes:

- `t`: time samples
- `x`: structured array of state variables
- `dx`: structured array of derivatives
- `aux`: structured array of auxiliary variables when present

A second `trajectory()` call with the default `update_x0=True` advances to the next window instead of rerunning the same one. See [continuation.md](continuation.md) when you want that handoff to be explicit.

## Ensemble inputs and layout

`set_ensemble(...)` accepts either:

- a mapping from variable or parameter names to scalars or one-dimensional arrays
- a full two-dimensional NumPy array with shape `(ensemble_size, num_variables)` or `(ensemble_size, num_parameters)`

clODE accepts these arrays directly, so most workflows do not need to think about storage details beyond the shapes above.

For more on RHS definitions, XPP conversion, auxiliary variables, and stochastic terms, see [specifying_odes.md](specifying_odes.md) and [api_reference.md](api_reference.md).
