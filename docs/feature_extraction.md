# Feature extraction

`FeatureSimulator` computes trajectory statistics on the device while integrating, which keeps memory use low even for large ensembles. This is the right tool when you want summary quantities such as periods, extrema, event counts, or event timestamps instead of full trajectories.

## Built-in observers

clODE currently ships the following observer modes through the public `Observer` enum:

- `Observer.summary`: summary statistics without event detection
- `Observer.local_max`: local maxima/minima events with IMI and amplitude readouts
- `Observer.threshold_crossing`: one-pass absolute threshold-crossing observer with runtime direction selection
- `Observer.normalized_threshold_crossing`: lean warmup-derived threshold-crossing observer with runtime direction selection
- `Observer.schmitt_trigger`: one-pass absolute Schmitt-trigger observer with up/down transition timestamps
- `Observer.normalized_schmitt_trigger`: warmup-normalized Schmitt-trigger observer with up/down transition timestamps
- `Observer.normalized_neighborhood_return`: warmup-normalized neighborhood-return observer with event timestamps and trajectory summaries

These are canonical names; there are no observer-name aliases in the enum.

The exact feature names depend on the observer. Use `get_feature_names()` on a configured simulator or `get_feature_names()` on the resulting `ObserverOutput` to inspect what is available.

## Example

The example below measures the period of the Van der Pol oscillator across an ensemble of `mu` values.

This example is also stored as [examples/van_der_pol_periods.py](https://github.com/patrickfletcher/clODE/blob/main/examples/van_der_pol_periods.py).

## Configuring an observer

Threshold-family observers now have a preferred semantic config surface through the `observer_configuration=` constructor argument and `set_observer_configuration(...)`.

Use `ThresholdCrossingConfig` when you want one threshold plus a direction selector. Pair it with `Observer.threshold_crossing` when that threshold should be interpreted in the units of `event_var`, or with `Observer.normalized_threshold_crossing` when the same `threshold` field should be interpreted as a warmup-derived fraction of the observed amplitude range:

```python
simulator = clode.FeatureSimulator(
    ...,
    observer=clode.Observer.normalized_threshold_crossing,
    observer_configuration=clode.ThresholdCrossingConfig(
        event_var="x",
        threshold=0.75,
        direction=clode.EventDirection.rising,
        min_amp=0.1,
        max_event_timestamps=16,
    ),
)
```

Use `SchmittTriggerConfig` when you want separate up/down boundaries on the Schmitt families. Pair it with `Observer.schmitt_trigger` for absolute thresholds or with `Observer.normalized_schmitt_trigger` for warmup-derived fractional thresholds:

```python
simulator.set_observer_configuration(
    clode.SchmittTriggerConfig(
        event_var="x",
        feature_var="x",
        x_up_threshold=0.3,
        x_down_threshold=0.2,
        min_amp=0.1,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.normalized_schmitt_trigger,
)
```

`x_up_threshold` must be greater than or equal to `x_down_threshold`. Equality is allowed and keeps the Schmitt-style up/down readouts while collapsing the hysteresis band to one shared boundary.

Use `LocalMaximumConfig` when you want local-maximum timestamps and values together with IMI, amplitude, and trajectory summaries:

```python
simulator.set_observer_configuration(
    clode.LocalMaximumConfig(
        event_var="x",
        max_event_timestamps=16,
    ),
    observer=clode.Observer.local_max,
)
```

Use `NeighborhoodReturnConfig` when you want a two-pass normalized neighborhood-return trigger:

```python
simulator.set_observer_configuration(
    clode.NeighborhoodReturnConfig(
        event_var="x",
        anchor_threshold=0.25,
        radius=0.15,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.normalized_neighborhood_return,
)
```

`anchor_threshold` chooses the falling `event_var` section that latches the sampled anchor point `x0`, and `radius` chooses how local the normalized full-state neighborhood stays around that anchor. This can be a good fit when a single scalar threshold is not distinctive enough to represent one recurrence cleanly. On a simple limit cycle it behaves like a light return-map trigger around one anchored point on the orbit; on a more complex or multi-lobed cycle it can still separate nearby passes because the exit rule is based on the full normalized state rather than only one scalar crossing.

When you call `set_observer_configuration(...)` with one of the shared config classes, pass `observer=` whenever you are switching families or when the simulator was not already constructed with the desired observer. The config type no longer implies absolute versus warmup-derived semantics on its own.

Single-observer configs such as `LocalMaximumConfig` and `NeighborhoodReturnConfig` are unambiguous, so clODE can infer their observer when you construct a simulator from the config alone. Passing `observer=` explicitly is still fine when you want the selection to stay obvious at the call site.

`set_observer_parameters(...)`, the public `ObserverParams` bundle, and the constructor `observer_*` keyword arguments remain available as compatibility surfaces when you need the older broad bundle or when a family-specific config object does not exist yet.

For summary-only workflows, `Observer.summary` also accepts an explicit `SummaryObserverSelection` at construction time or through `set_summary_selection(...)`.

Threshold and normalized-threshold observers always expose retained event times and event count; readout-selection surfaces are intentionally removed until a full cross-family observer UX is standardized.

Common compatibility-surface options include:

- `event_var`: which variable is used for event detection
- `feature_var`: which variable is used for feature readout
- `event_direction`: which crossing direction to accept in directional threshold-crossing observers such as `threshold_crossing` and `normalized_threshold_crossing`
- `max_event_count`: how many events to accumulate
- `max_event_timestamps`: how many event timestamps to retain
- `min_amp`, `min_imi`, `nhood_radius`, `x_up_threshold`, `x_down_threshold`, `dx_up_threshold`, `dx_down_threshold`, `eps_dx`

The current observer families intentionally cover different workflows:

- `Observer.local_max` uses `LocalMaximumConfig.event_var` (with `eVarIx == fVarIx`) and keeps IMI, amplitude, trajectory summaries, plus separate local-maximum and local-minimum event streams.
- `Observer.normalized_neighborhood_return` uses `NeighborhoodReturnConfig.event_var`, `NeighborhoodReturnConfig.anchor_threshold`, and `NeighborhoodReturnConfig.radius` to pin a warmup-scaled sampled full-state anchor and then store interpolated neighborhood-exit times plus event count.
- `Observer.threshold_crossing` uses `ThresholdCrossingConfig.threshold` as an absolute value in the units of `event_var`, and `ThresholdCrossingConfig.direction` chooses rising, falling, or either crossing through that one level. On the compatibility path, that same threshold still maps to `x_up_threshold`.
- `Observer.normalized_threshold_crossing` uses the same `ThresholdCrossingConfig` fields, but interprets `threshold` as a fraction of the warmup-pass amplitude range of `event_var`. It keeps the same lean single-stream timestamp-plus-count readout as `Observer.threshold_crossing` while adding warmup-derived scaling and `min_amp` gating.
- `Observer.schmitt_trigger` uses `SchmittTriggerConfig.x_up_threshold` and `SchmittTriggerConfig.x_down_threshold` as absolute thresholds in the units of `event_var`, and currently retains only up/down transition timestamps plus event count.
- `Observer.normalized_schmitt_trigger` uses the same `SchmittTriggerConfig` fields, but interprets its `x_*` thresholds as warmup-derived fractions of the observed event-variable range. It keeps the same lean up/down timestamp-plus-count readout as `Observer.schmitt_trigger` while adding warmup-derived scaling and `min_amp` gating.
- Keep threshold crossing and Schmitt triggering separate in the public mental model even though equal up/down thresholds can collapse to a degenerate single-boundary case internally. They answer different workflow questions and expose different readout bundles.
- `dx_up_threshold` and `dx_down_threshold` remain compatibility parameters in `ObserverParams`; semantic Schmitt observers do not use derivative gates.

The most directly useful current safeguards for oscillation-oriented observers are:

- `ThresholdCrossingConfig.direction` in `Observer.threshold_crossing` to select rising, falling, or either absolute crossing through one threshold value
- `ThresholdCrossingConfig.threshold` in `Observer.threshold_crossing` to set that absolute crossing value directly in state-variable units
- `ThresholdCrossingConfig.threshold` in `Observer.normalized_threshold_crossing` to reuse one threshold scalar while letting the live crossing level scale with the warmup-pass amplitude of `event_var`
- `ThresholdCrossingConfig.min_amp` and `ThresholdCrossingConfig.max_event_count` in either threshold-crossing family to suppress tiny oscillations and cap retained events without changing the observer schema
- `SchmittTriggerConfig.min_amp` and `SchmittTriggerConfig.max_event_count` in `Observer.schmitt_trigger` or `Observer.normalized_schmitt_trigger` to suppress small oscillations and bound event processing without changing the observer schema
- separate `SchmittTriggerConfig.x_up_threshold` and `SchmittTriggerConfig.x_down_threshold` values in the Schmitt families to set the hysteresis band in either absolute units or warmup-derived fractions, depending on the selected observer
- compatibility `ObserverParams` fields when you need broader low-level control for custom legacy-style tuning

Not every field is active in every built-in observer, so treat observer parameters as mode-specific rather than assuming every knob has the same effect everywhere. Recurring controls such as `min_amp` and `max_event_count` still live on the family configs where they actively shape oscillation-oriented outputs.

Persistent observer state also stays on the device for the duration of the solve, so its footprint depends on the observer mode, the model size, and `max_event_timestamps`. That internal state is distinct from the `ObserverOutput` readout object that `features()` returns. If you only need counts or summary statistics, keep `max_event_timestamps` as small as practical and prefer the lightest observer that answers the question.

For summary observers, the persistent state and output schema also scale with the selected summary groups. Changing a custom summary selection rebuilds the specialized OpenCL feature program, but it lets the stored summary state shrink to the requested subset instead of always following a one-variable or all-variable preset.

`Observer.schmitt_trigger` and `Observer.normalized_schmitt_trigger` store up/down transition times with inverse-linear interpolation of the active `x` boundary. The normalized family first converts its `x_*` thresholds from warmup-derived fractions of the observed `event_var` range into concrete live-pass thresholds; the absolute family uses configured values directly. `Observer.local_max` stores extrema using bounded three-sample quadratic refinement. `Observer.normalized_neighborhood_return` keeps sampled anchors and then refines each stored exit time by linearly interpolating the full normalized state between the last inside sample and the first outside sample of the neighborhood ball. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs.

`min_amp` is intentionally not identical across the one-pass absolute and two-pass normalized threshold families. In `Observer.threshold_crossing` and `Observer.schmitt_trigger` it is a live-pass range gate on `event_var`. In `Observer.normalized_threshold_crossing` and `Observer.normalized_schmitt_trigger` it is compared against the warmup-derived amplitude that defines the normalized thresholds. For steady-state workflows after a transient that difference is usually fine, but it is still part of the observer contract.

Compatibility-surface example for `Observer.threshold_crossing`:

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

`Observer.threshold_crossing` exposes one event stream, so `get_timestamps("event")` returns stored crossing times.

`Observer.local_max` keeps separate local-maximum and local-minimum streams together with IMI, amplitude, and trajectory-summary outputs. Use `get_timestamps("local maximum")` and `get_event_data("local maximum", type="value")` for maxima.

`Observer.normalized_threshold_crossing` also exposes one `event` stream, but it converts `ThresholdCrossingConfig.threshold` into one concrete boundary from the warmup-pass amplitude of `event_var` before the live pass starts.

`Observer.schmitt_trigger` and `Observer.normalized_schmitt_trigger` expose separate `up` and `down` streams plus an event count. The normalized family converts its configured thresholds into concrete event-variable boundaries during warmup; the absolute family uses the configured thresholds directly.

`Observer.normalized_neighborhood_return` exposes one `event` stream with interpolated normalized-ball exit times plus an event count.

Trajectory-summary feature names use model variable names in user-facing outputs (for example `max v`, `min v`, `max dv/dt`). `Observer.normalized_neighborhood_return` also exposes warmup/live geometry fields per state variable as `{var}0` and `range {var}` (for example `v0` and `range v`).

For trigger-geometry figures that mirror the live event rules directly, see `examples/visualize_events_threshold_crossing.py`, `examples/visualize_events_schmitt_trigger.py`, `examples/visualize_events_localmax.py`, and `examples/visualize_events_neighborhood_return.py`. The neighborhood-return example now matches the live sampled-anchor plus interpolated-exit rule exactly.

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
