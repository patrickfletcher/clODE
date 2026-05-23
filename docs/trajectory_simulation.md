# Trajectory simulation

`TrajectorySimulator` stores sampled trajectories for plotting, post-processing, and inspection.

## Example - FitzHugh-Nagumo oscillator

The following example simulates the FitzHugh-Nagumo oscillator with the fixed-step RK4 stepper.

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

    dV: float = V - V ** 3 / 3 - w + current
    dw: float = epsilon * (V + a - b * w)

    derivatives[0] = dV
    derivatives[1] = dw

a = 0.7
variables = {"V": 1.0, "w": 0.0}
parameters = {"a": a, "b": 0.8, "current": 0.0, "epsilon": 1.0 / 12.5}

simulator = clode.TrajectorySimulator(
    rhs_equation=fitzhugh_nagumo,
    variables=variables,
    parameters=parameters,
    stepper=clode.Stepper.rk4,
    t_span=(0, 200),
    dt=0.02,
)

ensemble_parameters = {"current": np.arange(0.0, 0.6, 0.1)}

simulator.set_ensemble(parameters=ensemble_parameters)

trajectories = simulator.trajectory()

plt.figure(figsize=(8, 6))
for index in range(len(trajectories)):
    label = f"I={ensemble_parameters['current'][index]:.1f}"
    plt.plot(trajectories[index].x["V"], trajectories[index].x["w"], label=label)
plt.xlabel("V")
plt.ylabel("w")
plt.legend()
plt.title("FitzHugh-Nagumo phase plane")
plt.show()

# Plot the time series
plt.figure(figsize=(8, 6))
for index in range(0, len(trajectories), 2):
    label = f"I={ensemble_parameters['current'][index]}"
    plt.plot(trajectories[index].t, trajectories[index].x["V"], label=label)
plt.xlabel("t")
plt.ylabel("V")
plt.legend()
plt.title("FitzHugh-Nagumo time series")
plt.show()
```

## Continuation and split windows

`trajectory()` continues the device state by default, but each returned `TrajectoryOutput`
contains only the samples from the current requested window.

To reproduce one long run with multiple windows:

- advance the next requested window from `get_final_time()` rather than the requested end time
- run `trajectory()` again
- concatenate the returned windows on the host and drop the duplicated boundary sample from the
    later window

See `continuation.md` for the full continuation semantics and `examples/continuation.py` for a
complete runnable example.
