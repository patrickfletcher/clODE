# Trajectory Simulation

`TrajectorySimulator` stores sampled trajectories for plotting, post-processing, and inspection.

## Understanding Device State and Outputs

When you call `trajectory()`, clODE returns only the trajectory samples generated during that call. Each call produces a new `TrajectoryOutput` object for that window only; it does not concatenate past windows for you.

With the default `update_x0=True`, repeated `trajectory()` calls continue both state and time into the next nominal window automatically. Pass `update_x0=False` when you want to inspect a window without applying that default continuation handoff afterward.

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

## Continuation and Split Windows

By default, repeated `trajectory()` calls already march forward window by window: after each solve, clODE promotes `xf -> x0`, continues the hidden timebase from the attained per-item final times, and shifts the nominal requested window by one duration.

To reproduce one long run using multiple windows:

1. Run `trajectory()` once with your first time window.
2. Run `trajectory()` again for each additional window.
3. Concatenate the returned windows on the host and drop the duplicate boundary sample from later windows if needed.

Use `set_tspan((start, end))` only when you intentionally want to reset or branch to a new shared absolute start time. `get_tspan()` reports the nominal requested window; `get_final_time()` reports the attained absolute final time for each ensemble member.

For a complete runnable example and detailed continuation semantics, see [continuation.md](continuation.md) and [examples/continuation.py](https://github.com/patrickfletcher/clODE/blob/main/examples/continuation.py).

## Solver diagnostics

After `trajectory()`, inspect solver diagnostics on the simulator instead of inferring them from the returned samples:

- `get_status()` reports whether the solve completed, hit `max_steps`, stopped at a terminal event, stopped because trajectory storage filled, or made no progress because the requested window collapsed in runtime precision
- `get_step_count()` reports accepted step counts
- `get_last_accepted_dt()` reports the width of the last accepted step
- `get_dt()` reports the continuation step size currently stored on the device, which for adaptive steppers is the next step size the controller would try on a continued solve rather than the last accepted width

This is especially useful when `TrajectorySimulator` stops early because of `max_store` or `nout`-driven storage policy rather than because the requested end time was reached.
