# Feature extraction

`FeatureSimulator` computes trajectory statistics on the device while integrating, which keeps memory use low even for large ensembles. This is the right tool when you want summary quantities such as periods, extrema, event counts, or event timestamps instead of full trajectories.

## Built-in observers

clODE currently ships the following observer modes through the public `Observer` enum:

- `Observer.summary`: selected summary statistics for explicit state, auxiliary, and slope subsets
- `Observer.basic`: basic summary statistics for one variable
- `Observer.basic_all_variables`: compatibility preset for the full summary bundle across all state and auxiliary variables
- `Observer.local_max`: local-maximum event tracking
- `Observer.threshold_1`: absolute threshold-crossing event tracking with runtime direction selection
- `Observer.neighbourhood_1`
- `Observer.neighbourhood_2`
- `Observer.threshold_2`: warmup-derived fractional Schmitt-trigger event tracking and period-style measurements

`Observer.basic` and `Observer.basic_all_variables` remain supported, but they now route through the same summary-observer family as `Observer.summary`.

The exact feature names depend on the observer. Use `get_feature_names()` on a configured simulator or `get_feature_names()` on the resulting `ObserverOutput` to inspect what is available.

## Example

The example below measures the period of the Van der Pol oscillator across an ensemble of `mu` values.

This example is also stored as [examples/van_der_pol_periods.py](https://github.com/patrickfletcher/clODE/blob/main/examples/van_der_pol_periods.py) so the docs and the runnable script stay aligned.

```py source run
--8<-- "examples/van_der_pol_periods.py"
```

## Configuring an observer

Most observer updates go through `set_observer_parameters(...)`. The public `ObserverParams` bundle remains available when you want to pass one object through the legacy-compatible surface.

For summary-only workflows, `Observer.summary` also accepts an explicit `SummaryObserverSelection` at construction time or through `set_summary_selection(...)`.

Common options include:

- `event_var`: which variable is used for event detection
- `feature_var`: which variable is used for feature readout
- `event_direction`: which crossing direction to accept in directional threshold-crossing observers such as `threshold_1`
- `max_event_count`: how many events to accumulate
- `max_event_timestamps`: how many event timestamps to retain
- `min_amp`, `min_imi`, `nhood_radius`, `x_up_threshold`, `x_down_threshold`, `dx_up_threshold`, `dx_down_threshold`, `eps_dx`

The current threshold observers intentionally cover different workflows:

- `threshold_1` uses `x_up_threshold` as an absolute value in the units of `event_var`, and `event_direction` chooses rising, falling, or either crossing through that one level.
- `threshold_2` interprets `x_up_threshold` and `x_down_threshold` as fractions of the warmup-pass amplitude range of `event_var`, so the effective thresholds scale with the observed oscillation range. That is often more robust across parameter sweeps where the absolute event-variable range changes.
- Treat `threshold_2` as the current Schmitt-style family: separate up/down fractions define the hysteresis band and support period, duty, and active-state measurements. Even though equal up/down fractions collapse to a degenerate single-boundary case internally, keeping threshold crossing and Schmitt triggering separate in the public mental model is clearer for configuration and readout selection.
- `dx_up_threshold` and `dx_down_threshold` are mainly useful as extra gates on noisy or stochastic traces. Smooth deterministic runs often do not need them.

The most directly useful current safeguards for oscillation-oriented observers are:

- `event_direction` in `threshold_1` to select rising, falling, or either absolute crossing through one threshold value
- `x_up_threshold` in `threshold_1` to set that absolute crossing value directly in state-variable units
- `min_amp` to suppress event measurement below a chosen amplitude floor
- separate `x_up_threshold` and `x_down_threshold` values in `threshold_2` to add warmup-derived Schmitt-trigger-style hysteresis and reduce chatter
- `dx_up_threshold` and `dx_down_threshold` in `threshold_2` when noisy shallow crossings need an additional slope gate

Not every field is active in every built-in observer, so treat observer parameters as mode-specific rather than assuming every knob has the same effect everywhere.

Persistent observer state also stays on the device for the duration of the solve, so its footprint depends on the observer mode, the model size, and `max_event_timestamps`. That internal state is distinct from the `ObserverOutput` readout object that `features()` returns. If you only need counts or summary statistics, keep `max_event_timestamps` as small as practical and prefer the lightest observer that answers the question.

For summary observers, the persistent state and output schema also scale with the selected summary groups. Changing a custom summary selection rebuilds the specialized OpenCL feature program, but it lets the stored summary state shrink to the requested subset instead of always following a one-variable or all-variable preset.

`threshold_2` converts its `x_up_threshold` and `x_down_threshold` inputs from warmup-pass fractions of the observed `event_var` amplitude into concrete boundaries, then stores up/down transition times with inverse-linear interpolation of the active boundary. When a `dx` threshold is zero, that slope gate is ignored; when it is nonzero, the stored transition time is the later of the active `x` and `dx` boundary crossings within the step. `local_max` stores extrema using bounded three-sample quadratic refinement. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

For `threshold_1`:

```python
integrator.set_observer_parameters(
    event_var="x",
    feature_var="x",
    event_direction=clode.EventDirection.rising,
    x_up_threshold=-40.0,
    max_event_timestamps=16,
)
```

Custom summary-selection example:

```python
selection = clode.SummaryObserverSelection(
    state={"y": ("max", "mean")},
    slope={"x": "min"},
    aux={"sum": "mean"},
)

simulator = clode.FeatureSimulator(
    ...,
    observer=clode.Observer.summary,
    summary_selection=selection,
)
```

## Reading results

`features()` returns an `ObserverOutput` readout object. Common accessors are:

- `get_feature_names()`
- `get_var_mean(name)`
- `get_var_min(name)`
- `get_var_max(name)`
- `get_var_mean_slope(name)`, `get_var_min_slope(name)`, `get_var_max_slope(name)`
- `get_var_count(name)`
- `get_event_data(name, type="time")`

`Observer.threshold_1` exposes one `threshold` event stream, so `get_timestamps("threshold")` returns the stored absolute crossing times for that observer.

`Observer.threshold_2` exposes separate `up` and `down` streams after converting its warmup-derived fractional thresholds into concrete event-variable boundaries for the live pass.

Solver diagnostics stay on the simulator rather than in `ObserverOutput`:

- `get_status()` returns solver-owned completion or early-stop codes
- `get_step_count()` returns accepted step counts
- `get_last_accepted_dt()` returns the width of the last accepted step
- `get_dt()` returns the continuation step size stored on the device; for adaptive steppers that can differ from `get_last_accepted_dt()` because it is the next step size the controller would attempt on a continued solve

If you need stop reasons or stepping diagnostics, query the simulator after `features()` rather than expecting observer feature names such as step-count or `dt` summaries.

## Continuation and repeated calls

`features()` continues device state by default and also continues the observer state unless
it is explicitly reinitialized. The requested `t_span` does not advance automatically.

For exact absolute-time continuation, advance the next requested window from
`get_final_time()` before calling `features()` again. This is especially important for
fixed-step methods, time-based feature accumulators, event timestamps, and non-autonomous
systems.

For autonomous systems, prefer feature windows whose local `t_span` starts near `0` when
absolute time is not part of the model. Large absolute times still coarsen stored float32
absolute timestamps even though the observers keep elapsed-time statistics separate from
those large absolute values.

See `continuation.md` for the full continuation model and `examples/continuation.py` for a
runnable comparison between one long feature run and split-window continuation.

For empirical float32 demonstrations of amplitude floors, threshold hysteresis, derivative thresholds, compensated means, elapsed-time origins, and time accumulation limits, see [numerical_accuracy.md](numerical_accuracy.md).

## Custom observers

The built-in observers are stable and supported. Custom observer authoring is not part of the public API.
