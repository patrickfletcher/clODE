# Feature extraction

`FeatureSimulator` computes trajectory statistics on the device while integrating, which keeps memory use low even for large ensembles. This is the right tool when you want summary quantities such as periods, extrema, event counts, or event timestamps instead of full trajectories.

## Built-in observers

clODE currently ships the following observer modes through the public `Observer` enum:

- `Observer.summary`: selected summary statistics for explicit state, auxiliary, and slope subsets
- `Observer.basic`: basic summary statistics for one variable
- `Observer.basic_all_variables`: compatibility preset for the full summary bundle across all state and auxiliary variables
- `Observer.local_extremum`: lean local-extremum event tracking with selectable polarity
- `Observer.local_max`: retained legacy fully featured maxima-oriented observer
- `Observer.threshold_crossing`: one-pass absolute threshold-crossing observer with runtime direction selection
- `Observer.normalized_threshold_crossing`: lean warmup-derived threshold-crossing observer with runtime direction selection
- `Observer.schmitt_trigger`: one-pass absolute Schmitt-trigger observer with retained up/down timestamps
- `Observer.normalized_schmitt_trigger`: lean warmup-derived Schmitt-trigger observer with retained up/down timestamps
- `Observer.threshold_2`: retained legacy fully featured normalized Schmitt-trigger observer
- `Observer.neighborhood_return`: lean warmup-derived normalized neighborhood-return observer
- `Observer.neighbourhood_1`
- `Observer.neighbourhood_2`: retained legacy fully featured normalized neighborhood-return observer

`Observer.basic` and `Observer.basic_all_variables` remain supported, but they now route through the same summary-observer family as `Observer.summary`.

`Observer.local_max`, `Observer.threshold_2`, and `Observer.neighbourhood_2` remain supported as retained legacy fully featured workflows.

The exact feature names depend on the observer. Use `get_feature_names()` on a configured simulator or `get_feature_names()` on the resulting `ObserverOutput` to inspect what is available.

## Example

The example below measures the period of the Van der Pol oscillator across an ensemble of `mu` values.

This example is also stored as [examples/van_der_pol_periods.py](https://github.com/patrickfletcher/clODE/blob/main/examples/van_der_pol_periods.py) so the docs and the runnable script stay aligned.

```py source run
--8<-- "examples/van_der_pol_periods.py"
```

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

Use `SchmittTriggerConfig` when you want separate up/down boundaries and optional slope gates. Pair it with `Observer.schmitt_trigger` for absolute thresholds or with `Observer.normalized_schmitt_trigger` for warmup-derived fractional thresholds. Use `Observer.threshold_2` when you need the retained legacy fully featured normalized Schmitt readout:

```python
simulator.set_observer_configuration(
    clode.SchmittTriggerConfig(
        event_var="x",
        feature_var="x",
        x_up_threshold=0.3,
        x_down_threshold=0.2,
        dx_up_threshold=0.0,
        dx_down_threshold=0.0,
        min_amp=0.1,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.normalized_schmitt_trigger,
)
```

Use `LocalExtremumConfig` when you want lean local-extremum timestamps and values without the broader `Observer.local_max` summary bundle:

```python
simulator.set_observer_configuration(
    clode.LocalExtremumConfig(
        variable="x",
        polarity=clode.ExtremumPolarity.maximum,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.local_extremum,
)
```

Use `NeighborhoodReturnConfig` when you want the lean two-pass normalized neighborhood-return trigger without the broader `Observer.neighbourhood_2` readout bundle:

```python
simulator.set_observer_configuration(
    clode.NeighborhoodReturnConfig(
        event_var="x",
        anchor_threshold=0.25,
        radius=0.15,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.neighborhood_return,
)
```

When you call `set_observer_configuration(...)` with one of the shared config classes, pass `observer=` whenever you are switching families or when the simulator was not already constructed with the desired observer. The config type no longer implies absolute versus warmup-derived semantics on its own.

Single-observer configs such as `LocalExtremumConfig` and `NeighborhoodReturnConfig` are unambiguous, so clODE can infer their observer when you construct a simulator from the config alone. Passing `observer=` explicitly is still fine when you want the selection to stay obvious at the call site.

`set_observer_parameters(...)`, the public `ObserverParams` bundle, and the constructor `observer_*` keyword arguments remain available as compatibility surfaces when you need the older broad bundle or when a family-specific config object does not exist yet.

For summary-only workflows, `Observer.summary` also accepts an explicit `SummaryObserverSelection` at construction time or through `set_summary_selection(...)`.

Common compatibility-surface options include:

- `event_var`: which variable is used for event detection
- `feature_var`: which variable is used for feature readout
- `event_direction`: which crossing direction to accept in directional threshold-crossing observers such as `threshold_crossing` and `normalized_threshold_crossing`
- `max_event_count`: how many events to accumulate
- `max_event_timestamps`: how many event timestamps to retain
- `min_amp`, `min_imi`, `nhood_radius`, `x_up_threshold`, `x_down_threshold`, `dx_up_threshold`, `dx_down_threshold`, `eps_dx`

The current observer families intentionally cover different workflows:

- `Observer.local_extremum` uses `LocalExtremumConfig.variable` and `LocalExtremumConfig.polarity` to store a lean stream of local-extremum timestamps and event values plus event count.
- `Observer.local_max` remains the legacy fully featured maxima-oriented workflow. It keeps the broader IMI, amplitude, and trajectory-summary bundle together with separate `localmax` and `localmin` event streams.
- `Observer.neighborhood_return` uses `NeighborhoodReturnConfig.event_var`, `NeighborhoodReturnConfig.anchor_threshold`, and `NeighborhoodReturnConfig.radius` to pin a warmup-scaled full-state anchor and then store neighborhood-exit times plus event count.
- `Observer.neighbourhood_2` remains the legacy fully featured normalized neighborhood-return workflow. It keeps the broader period, peaks, and trajectory-summary bundle around the same general trigger idea.

- `Observer.threshold_crossing` uses `ThresholdCrossingConfig.threshold` as an absolute value in the units of `event_var`, and `ThresholdCrossingConfig.direction` chooses rising, falling, or either crossing through that one level. On the compatibility path, that same threshold still maps to `x_up_threshold`.
- `Observer.normalized_threshold_crossing` uses the same `ThresholdCrossingConfig` fields, but interprets `threshold` as a fraction of the warmup-pass amplitude range of `event_var`. It keeps the same lean single-stream timestamp-plus-count readout as `Observer.threshold_crossing` while adding warmup-derived scaling and `min_amp` gating.
- `Observer.schmitt_trigger` uses `SchmittTriggerConfig.x_up_threshold` and `SchmittTriggerConfig.x_down_threshold` as absolute thresholds in the units of `event_var`, and currently retains only up/down transition timestamps plus event count.
- `Observer.normalized_schmitt_trigger` uses the same `SchmittTriggerConfig` fields, but interprets its `x_*` and optional `dx_*` thresholds as warmup-derived fractions of the observed event-variable range. It keeps the same lean up/down timestamp-plus-count readout as `Observer.schmitt_trigger` while adding warmup-derived scaling and `min_amp` gating.
- `Observer.threshold_2` remains the legacy fully featured normalized Schmitt workflow. It uses the same warmup-derived threshold semantics as `Observer.normalized_schmitt_trigger`, but it also retains the broader period, duty, active-dip, and trajectory-summary bundle.
- Keep threshold crossing and Schmitt triggering separate in the public mental model even though equal up/down thresholds can collapse to a degenerate single-boundary case internally. They answer different workflow questions and expose different readout bundles.
- `dx_up_threshold` and `dx_down_threshold` are mainly useful as extra gates on noisy or stochastic traces. Smooth deterministic runs often do not need them.

The most directly useful current safeguards for oscillation-oriented observers are:

- `ThresholdCrossingConfig.direction` in `Observer.threshold_crossing` to select rising, falling, or either absolute crossing through one threshold value
- `ThresholdCrossingConfig.threshold` in `Observer.threshold_crossing` to set that absolute crossing value directly in state-variable units
- `ThresholdCrossingConfig.threshold` in `Observer.normalized_threshold_crossing` to reuse one threshold scalar while letting the live crossing level scale with the warmup-pass amplitude of `event_var`
- `ThresholdCrossingConfig.min_amp` and `ThresholdCrossingConfig.max_event_count` in either threshold-crossing family to suppress tiny oscillations and cap retained events without changing the observer schema
- `SchmittTriggerConfig.min_amp` and `SchmittTriggerConfig.max_event_count` in `Observer.schmitt_trigger`, `Observer.normalized_schmitt_trigger`, or `Observer.threshold_2` to suppress small oscillations and bound event processing without changing the observer schema
- separate `SchmittTriggerConfig.x_up_threshold` and `SchmittTriggerConfig.x_down_threshold` values in the Schmitt families to set the hysteresis band in either absolute units or warmup-derived fractions, depending on the selected observer
- `SchmittTriggerConfig.dx_up_threshold` and `SchmittTriggerConfig.dx_down_threshold` in the Schmitt families when noisy shallow crossings need an additional slope gate; `dx_down_threshold` is specified as a positive magnitude for the required negative downward slope

Not every field is active in every built-in observer, so treat observer parameters as mode-specific rather than assuming every knob has the same effect everywhere. Recurring controls such as `min_amp` and `max_event_count` currently live on the family configs where they actively shape oscillation-oriented outputs; clODE does not yet add a separate shared oscillation bundle for them.

Persistent observer state also stays on the device for the duration of the solve, so its footprint depends on the observer mode, the model size, and `max_event_timestamps`. That internal state is distinct from the `ObserverOutput` readout object that `features()` returns. If you only need counts or summary statistics, keep `max_event_timestamps` as small as practical and prefer the lightest observer that answers the question.

For summary observers, the persistent state and output schema also scale with the selected summary groups. Changing a custom summary selection rebuilds the specialized OpenCL feature program, but it lets the stored summary state shrink to the requested subset instead of always following a one-variable or all-variable preset.

`Observer.schmitt_trigger`, `Observer.normalized_schmitt_trigger`, and `Observer.threshold_2` all store up/down transition times with inverse-linear interpolation of the active boundary. The normalized semantic family and `Observer.threshold_2` first convert their `x_*` and `dx_*` fields from warmup-derived fractions of the observed `event_var` range into concrete live-pass thresholds; the absolute family uses the configured values directly. When a `dx` threshold is zero, that slope gate is ignored; when it is nonzero, the stored transition time is the later of the active `x` and `dx` boundary crossings within the step. `Observer.threshold_2` additionally keeps the heavier legacy period, duty, and active-dip bundle. `Observer.local_extremum` and `Observer.local_max` store extrema using bounded three-sample quadratic refinement. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

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

`Observer.threshold_crossing` exposes one `threshold` event stream, so `get_timestamps("threshold")` returns the stored absolute crossing times for that observer.

`Observer.local_extremum` exposes one `local_extremum` event stream with timestamps and event values plus an event count. Use `get_timestamps("local_extremum")` for times and `get_event_data("local_extremum", type="value")` for the stored values.

`Observer.local_max` remains the heavier maxima-oriented legacy workflow. It keeps separate `localmax` and `localmin` event streams together with IMI, amplitude, and trajectory-summary outputs.

`Observer.normalized_threshold_crossing` also exposes one `threshold` event stream, but it converts `ThresholdCrossingConfig.threshold` into one concrete boundary from the warmup-pass amplitude of `event_var` before the live pass starts.

`Observer.schmitt_trigger` and `Observer.normalized_schmitt_trigger` expose separate `up` and `down` streams plus an event count. The normalized family converts its configured thresholds into concrete event-variable boundaries during warmup; the absolute family uses the configured thresholds directly.

`Observer.threshold_2` exposes the same `up` and `down` streams while also retaining the legacy heavier period, duty, active-dip, and trajectory-summary readout bundle.

`Observer.neighborhood_return` exposes one `neighborhood_return` event stream with normalized neighborhood-exit times plus an event count.

`Observer.neighbourhood_2` keeps the same broad trigger family but also retains the legacy heavier period, peaks, and trajectory-summary readout bundle.

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
