# Continuation and repeated runs

Repeated calls on the same simulator are stateful. With the default `update_x0=True`, clODE treats each new solve as the next chunk of the same run: it promotes the previous final state to the next initial state, advances the hidden solver timebase from the attained final times, and shifts the public `t_span` forward by one nominal duration.

## What continues automatically

By default, repeated calls preserve and advance solver-side state:

- `transient()`, `trajectory()`, and `features()` all continue both state and time when `update_x0=True`
- device `dt` is preserved as continuation or controller state
- RNG state is preserved
- `FeatureSimulator` also preserves persistent observer state unless you explicitly reinitialize it

If you pass `update_x0=False`, clODE still runs the requested window, but it does not perform that default state-and-time handoff afterward.

## Requested window versus attained time

Two public time views are intentionally different:

- `get_tspan()` returns the nominal requested window stored on the simulator
- `get_final_time()` returns the attained absolute final time for each ensemble member from the last solve

They can differ after fixed-step overshoot, early termination, or diverged ensemble timing. In normal continued solves, clODE keeps using the attained per-item times internally even though `get_tspan()` still reports only the nominal window you requested.

## The explicit controls

- `set_tspan((start, end))` resets both the nominal requested window and the hidden solver timebase to `start`
- `shift_x0()` promotes `xf -> x0` without changing the nominal requested window
- `shift_tspan()` advances the nominal requested window by one duration and continues the solver-owned timebase from the attained per-item `tf`

These helpers are useful when you want to split the state and time handoff yourself instead of relying on the default repeated-call behavior.

## Common workflows

### Keep solving forward in equal windows

For the common case, just call the same solve method again:

```python
simulator.features()
simulator.features()
simulator.features()
```

With the default `update_x0=True`, each call continues the previous one.

### Reset or branch to a new shared start time

Use `set_tspan()` when you want the next solve to use a specific shared absolute start time:

```python
simulator.set_tspan((200.0, 250.0))
simulator.features(update_x0=False)
```

This changes the timebase only. Use it together with a fresh simulator or new initial conditions when you also want a full state reset.

### Manage the handoff explicitly

If you need host-side logic between solves, disable the automatic handoff and apply the pieces yourself:

```python
simulator.transient(update_x0=False)
simulator.shift_x0()
simulator.shift_tspan()
```

Apply only one of those helpers when that is the behavior you want.

## Trajectory windows

`TrajectorySimulator.trajectory()` returns only the samples from the current window. If you split a long run into multiple windows, concatenate the returned `TrajectoryOutput` objects on the host and drop the duplicated boundary sample from later windows when needed.

## Feature windows and observer state

Repeated `features()` calls continue persistent observer state by default, which is what you want when multiple windows should behave like one long observation window.

When you want a fresh observer pass on the current state instead, rerun the observer initialization path with `initialize_observer=True` or start from a fresh simulator configuration.

## Solver diagnostics after a run

- `get_status()` returns the solver-owned stop reason
- `get_step_count()` returns accepted step counts
- `get_last_accepted_dt()` returns the width of the last accepted step
- `get_dt()` returns the continuation step size currently stored on the device

For adaptive steppers, `get_dt()` and `get_last_accepted_dt()` can differ because the controller proposes the next step size after each accepted step.

## Single-precision note

clODE keeps solver time in compensated solve-relative form and lets observers consume that solver-owned elapsed time instead of rebuilding large float32 time differences on every step. That materially improves endpoint accuracy, periods, durations, and time-weighted summaries in long single-precision runs, especially when the solve window starts at a large absolute time.

It does not change float32 spacing itself. Absolute timestamps can still quantize in `ulp(t)`-sized jumps, and double precision or shorter windows remain the safer choice when fine absolute-time fidelity is the requirement.

For the detailed demonstrations and tradeoffs, see [numerical_accuracy.md](numerical_accuracy.md).

## Runnable example

The repository includes a full example script in `examples/continuation.py` that compares a single long run against split-window continuation for both `TrajectorySimulator` and `FeatureSimulator`.
