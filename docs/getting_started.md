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

Repeated runs are stateful by default, but exact continuation also depends on how
the requested `t_span` is advanced between calls. See `continuation.md` before
building long-running split-window workflows.

## Example: Van der Pol period measurement

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

The feature example below measures oscillation period for an ensemble of `mu` values without storing full trajectories.

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
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=t_span,
)

integrator.set_ensemble(parameters={"mu": mu_values})
integrator.transient()

observer_output = integrator.features()
print(observer_output.get_var_mean("period"))
```

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
    stepper=clode.Stepper.dormand_prince,
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

## Ensemble inputs and layout

`set_ensemble(...)` accepts either:

- a mapping from variable or parameter names to scalars or one-dimensional arrays
- a full two-dimensional NumPy array with shape `(ensemble_size, num_variables)` or `(ensemble_size, num_parameters)`

Internally, clODE flattens problem data in column-major order before sending it to OpenCL buffers. That layout is part of the backend implementation and normally does not need to be handled directly in user code.

For more on RHS definitions, XPP conversion, auxiliary variables, and stochastic terms, see [specifying_odes.md](specifying_odes.md) and [api_reference.md](api_reference.md).
