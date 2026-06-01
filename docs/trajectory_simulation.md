# Trajectory Simulation

`TrajectorySimulator` stores sampled trajectories for plotting, post-processing, and inspection.

## Understanding Device State and Outputs

When you call `trajectory()`, clODE returns only the trajectory samples generated during that call. However, the **solver state persists on the device** between calls. This is important for continuation: if you run `trajectory()` a second time without explicitly advancing the time window, the solver will resume from where it left off, not restart.

Each call to `trajectory()` returns a new `TrajectoryOutput` object containing only the samples from that particular window—not the samples from all previous calls. The device tracks where the solver is in continuous time, so you need to explicitly advance the requested time window before calling `trajectory()` again if you want to continue from the attained final time.

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

By default, `trajectory()` continues the device state and returns only the samples from the current requested window. To reproduce one long run using multiple windows:

1. Run `trajectory()` once with your first time window.
2. Get the attained final time from the device: `final_time = simulator.get_final_time()`.
3. For the next window, set the requested `t_span` to start from that attained final time: `simulator.set_tspan((final_time, final_time + window_length))`.
4. Run `trajectory()` again. The solver continues from `final_time` and returns only the new samples.
5. Concatenate the returned windows and drop the duplicate boundary sample from later windows if needed.

This is more accurate than using the originally requested end time, especially for adaptive steppers where the attained final time may differ from what was requested.

For a complete runnable example and detailed continuation semantics, see [continuation.md](continuation.md) and [examples/continuation.py](https://github.com/patrickfletcher/clODE/blob/main/examples/continuation.py).

## Solver diagnostics

After `trajectory()`, inspect solver diagnostics on the simulator instead of inferring them from the returned samples:

- `get_status()` reports whether the solve completed, hit `max_steps`, stopped at a terminal event, stopped because trajectory storage filled, or made no progress because the requested window collapsed in runtime precision
- `get_step_count()` reports accepted step counts
- `get_last_accepted_dt()` reports the width of the last accepted step
- `get_dt()` reports the continuation step size currently stored on the device, which for adaptive steppers is the next step size the controller would try on a continued solve rather than the last accepted width

This is especially useful when `TrajectorySimulator` stops early because of `max_store` or `nout`-driven storage policy rather than because the requested end time was reached.
