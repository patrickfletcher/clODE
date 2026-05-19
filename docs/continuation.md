# Continuation and repeated runs

clODE supports repeated calls on the same simulator object. For users, the important
question is what state persists across calls and what has to be advanced explicitly.

## What persists across calls

By default, repeated calls preserve solver-side state:

- `x0` is replaced with the previous `xf` when `update_x0=True`
- device `dt` is preserved
- RNG state is preserved
- `FeatureSimulator` also preserves observer state unless the observer is reinitialized

What does **not** change automatically is the requested time window. The next kernel
still starts from whatever `t_span` is currently configured on the simulator.

## Common continuation patterns

There are two common workflows:

- requested-window continuation: continue the solver state while reusing or shifting the
  nominal window you asked for
- exact absolute-time continuation: make the next window start from the attained final
  time of the previous run

For autonomous systems, state continuation is often the main thing you care about.
For non-autonomous systems, exact absolute time matters because the RHS depends on `t`.

## Fixed-step overshoot and `get_final_time()`

Adaptive steppers such as Dormand-Prince target the requested end time exactly. Fixed-step
methods do not: the attained final time can overshoot `t_span[1]` by up to one step.

For exact continuation across split windows, use the attained final time from
`get_final_time()` rather than the requested endpoint. This matters for:

- non-autonomous RHS evaluations
- observer timestamps and elapsed-time statistics
- exact parity between one long run and multiple split runs

The convenience method `shift_tspan()` advances the requested window by its requested
duration. That is still useful for requested-window continuation, but it is not guaranteed
to produce exact absolute-time continuation after a fixed-step run or an early stop.

## What to do in practice

For exact split-window continuation:

1. run `transient()`, `features()`, or `trajectory()`
2. read the attained end time from `get_final_time()`
3. set the next requested window from that attained time
4. run the next window

## Trajectory continuation

`TrajectorySimulator.trajectory()` returns only the samples from the current requested
window. If you split a long run into windows, concatenate the returned `TrajectoryOutput`
objects on the host and drop the duplicated boundary sample from the later window.

## Feature continuation

Repeated `features()` calls continue the observer state by default. Exact continuation of
time-based feature accumulators and event timestamps therefore requires the next requested
window to start from the previous attained final time.

If you call `features()` repeatedly without advancing `t_span`, you are not asking for the
same thing as a single long run. For time-based observers, that can make the accumulated
statistics inconsistent.

## A user-land helper pattern

For the common case where the ensemble shares one attained final time, the following helper
keeps the next requested window aligned with that attained time:

clODE intentionally leaves this as user code for now. Once ensemble members finish at
different times, there is no single correct shared-window update for the library to apply
automatically.

```python
import numpy as np


def advance_window_to_attained_final_time(simulator) -> None:
    start, end = simulator.get_tspan()
    duration = end - start
    final_times = np.asarray(simulator.get_final_time(), dtype=np.float64).reshape(-1)

    if not np.allclose(final_times, final_times[0], atol=1e-12, rtol=0.0):
        raise ValueError(
            "Shared t_span continuation needs a single agreed final time across the ensemble"
        )

    next_start = float(final_times[0])
    simulator.set_tspan((next_start, next_start + duration))
```

If ensemble members finish at different times, the caller has to choose an explicit policy.
There is no single correct shared `t_span` update in that case, so clODE does not yet ship
a built-in helper for it.

## Current limitation

clODE still advances fixed-step time with `ti += dt` inside the kernels. That means very long
absolute-time runs with very small `dt` eventually hit floating-point resolution limits. Once
`dt` is smaller than the spacing between adjacent representable values near `ti`, time updates
lose significance and the simulation can break down.

Implications:

- for autonomous systems, splitting a long run into shorter windows can delay the problem
- for non-autonomous systems, this remains a real numerical limitation

## Runnable example

The repository includes a full example script in `examples/continuation.py` that compares a
single long run against split-window continuation for both `TrajectorySimulator` and
`FeatureSimulator`.
