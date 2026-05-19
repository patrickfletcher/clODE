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
2. call `advance_tspan_to_attained_final_time()` when the ensemble shares one attained final time
3. otherwise, read `get_final_time()` and choose an explicit next-window policy in user code
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

## Built-in exact-continuation helper

For the common case where the ensemble shares one attained final time, use
`advance_tspan_to_attained_final_time()`:

```python
simulator.transient()
simulator.advance_tspan_to_attained_final_time()
simulator.transient()
```

This keeps the current requested duration but moves the next window start to the attained
`tf` from the previous solve. Unlike `shift_tspan()`, it uses the attained final time rather
than the requested endpoint.

If ensemble members finish at different times, `advance_tspan_to_attained_final_time()`
raises `ValueError`. In that case the caller still has to choose an explicit policy with
`get_final_time()` and `set_tspan()` because there is no single correct shared-window update.

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
