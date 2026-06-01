# Observers

Observers are the mechanism behind `FeatureSimulator`. They maintain per-ensemble state on the device and reduce trajectory information as integration proceeds, which lets clODE compute summary quantities without storing the full trajectory.

## Built-in observer modes

The current public observer modes are:

- `clode.Observer.summary`
- `clode.Observer.local_max`
- `clode.Observer.threshold_crossing`
- `clode.Observer.normalized_threshold_crossing`
- `clode.Observer.schmitt_trigger`
- `clode.Observer.normalized_schmitt_trigger`
- `clode.Observer.normalized_neighborhood_return`

These are canonical enum names; there are no observer-name aliases.

Choose an observer when constructing a `FeatureSimulator`:

```python
import clode


integrator = clode.FeatureSimulator(
    src_file="test/van_der_pol_oscillator.cl",
    variables={"x": 0.0, "y": 1.0},
    parameters={"mu": 1.0},
    observer=clode.Observer.schmitt_trigger,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
)
```

## Choosing a threshold workflow

The threshold catalog is easiest to reason about as two independent axes: one threshold versus Schmitt hysteresis, and absolute units versus warmup-derived fractions.

- `clode.Observer.threshold_crossing`: one-pass directional threshold crossing in the units of `event_var`
- `clode.Observer.normalized_threshold_crossing`: two-pass directional threshold crossing where the same `threshold` field is interpreted as a warmup-derived fraction of the observed amplitude range
- `clode.Observer.schmitt_trigger`: one-pass absolute Schmitt trigger with separate up/down boundaries in the units of `event_var`
- `clode.Observer.normalized_schmitt_trigger`: two-pass Schmitt trigger where `x_up_threshold` and `x_down_threshold` are interpreted as warmup-derived fractions in `[0, 1]`

Threshold crossing and Schmitt triggering stay separate user-facing concepts even though equal up/down thresholds can collapse to a degenerate single-boundary case internally. They expose different readout bundles and answer different workflow questions.

## Choosing an extremum or neighborhood workflow

- `clode.Observer.local_max`: canonical maxima-oriented workflow with IMI, amplitude, all-state summaries, and separate local-maximum/local-minimum event streams.
- `clode.Observer.normalized_neighborhood_return`: lean two-pass normalized neighborhood-return detector. Warmup fixes the normalization range, the live pass anchors on the first sampled point where `event_var` drops below `anchor_threshold`, and the observer stores linearly refined neighborhood-exit times plus event count, period/maxima/amplitude aggregates, and per-variable anchor/range fields.

## Observer parameters

Threshold-family observers now have preferred semantic config objects. Use `observer_configuration=` at construction time or `set_observer_configuration(...)` when you want a family-specific knob set instead of the broad compatibility bundle.

Use `clode.ThresholdCrossingConfig` for both threshold-crossing families:

```python
integrator = clode.FeatureSimulator(
    ...,
    observer=clode.Observer.normalized_threshold_crossing,
    observer_configuration=clode.ThresholdCrossingConfig(
        event_var="x",
        feature_var="x",
        threshold=0.75,
        direction=clode.EventDirection.rising,
        min_amp=0.1,
        max_event_timestamps=16,
    ),
)
```

Across the event-triggering semantic families, `event_var` selects the trigger geometry and the `min_amp` gate. `feature_var` selects the extrema/amplitude channel on families that split trigger geometry from measurement, which keeps threshold, Schmitt, and neighborhood-return aligned when you want one state variable to trigger events and another to supply oscillation readouts.

Use `clode.SchmittTriggerConfig` for the semantic Schmitt families:

```python
integrator.set_observer_configuration(
    clode.SchmittTriggerConfig(
        event_var="x",
        feature_var="x",
        max_event_count=128,
        max_event_timestamps=16,
        min_amp=0.1,
        x_up_threshold=0.3,
        x_down_threshold=0.2,
    ),
    observer=clode.Observer.normalized_schmitt_trigger,
)
```

`x_up_threshold` must be greater than or equal to `x_down_threshold`. Equality is allowed and collapses the hysteresis band to one shared boundary while keeping the Schmitt-style up/down readout streams.

Use `clode.LocalMaximumConfig` with `clode.Observer.local_max` for maxima-triggered extrema tracking:

```python
integrator.set_observer_configuration(
    clode.LocalMaximumConfig(
        event_var="x",
        min_amp=0.1,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.local_max,
)
```

`LocalMaximumConfig` stays one-channel by design: the same resolved variable drives extrema detection, stored extrema values, and the `min_amp` gate.

Use `clode.NeighborhoodReturnConfig` with `clode.Observer.normalized_neighborhood_return` when you want the lean normalized neighborhood-return trigger:

```python
integrator.set_observer_configuration(
    clode.NeighborhoodReturnConfig(
        event_var="x",
        feature_var="x",
        anchor_threshold=0.25,
        radius=0.15,
        min_amp=0.1,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.normalized_neighborhood_return,
)
```

`anchor_threshold` picks the falling `event_var` section that latches the anchor point `x0`, while `radius` sets the size of the normalized full-state neighborhood around that anchor. `event_var` also owns the family `min_amp` gate. `feature_var` selects the extrema/amplitude channel for the neighborhood observer; leave it equal to `event_var` for one-channel behavior, or point it at another state variable when the return trigger and oscillation readouts should differ. This can be a good fit when a one-variable threshold crossing is too ambiguous to identify one recurrence cleanly. On a simple limit cycle it behaves like a light-weight return-map trigger around one anchored point on the orbit; on a more complex or multi-lobed cycle it can still distinguish nearby passes because the exit test uses the full normalized state rather than only the scalar `event_var`.

When you call `set_observer_configuration(...)` with one of the shared config classes, pass `observer=` whenever you are switching families or when the simulator was not already constructed with the desired observer. The config class alone does not determine absolute versus warmup-derived interpretation.

Single-observer configs such as `clode.LocalMaximumConfig` and `clode.NeighborhoodReturnConfig` are unambiguous, so clODE can infer their observer when you construct a simulator from the config alone. Passing `observer=` explicitly is still useful when you want the selection to be obvious at the call site.

`clode.Observer.threshold_crossing` and `clode.Observer.normalized_threshold_crossing` always expose retained event times and event count. Readout-selection surfaces are intentionally removed until a broader cross-family observer UX is standardized.

Tune observer behavior with `set_observer_parameters(...)` when you want the older compatibility surface or need to pass the broad `ObserverParams` bundle.

```python
integrator.set_observer_parameters(
    event_var="x",
    feature_var="x",
    max_event_count=128,
    max_event_timestamps=16,
    min_amp=0.1,
    x_up_threshold=0.3,
    x_down_threshold=0.2,
)
```

Changing `max_event_timestamps` changes observer storage requirements and may trigger a rebuild of the OpenCL program.

Not every observer field affects every built-in observer. For the current event-triggering families, the most directly useful safeguards are:

- `min_amp` for suppressing oscillation/event measurements below a chosen amplitude floor
- `max_event_count` for terminating the feature solve after a bounded number of accepted events
- `threshold` plus `direction` in `clode.ThresholdCrossingConfig` for the threshold-crossing families, with the selected observer deciding whether that threshold is absolute or warmup-derived
- separate `x_up_threshold` and `x_down_threshold` values in `clode.SchmittTriggerConfig` for the Schmitt families, with the selected observer deciding whether those values are absolute or warmup-derived fractions
- `anchor_threshold` plus `radius` in `clode.NeighborhoodReturnConfig` for sampled-anchor neighborhood-return workflows
- semantic Schmitt-trigger families (`clode.Observer.schmitt_trigger` and `clode.Observer.normalized_schmitt_trigger`) when you need explicit up/down hysteresis around a threshold band

See [numerical_accuracy.md](numerical_accuracy.md) for empirical examples and practical tuning guidance.

`clode.Observer.schmitt_trigger` and `clode.Observer.normalized_schmitt_trigger` store up/down transition times with inverse-linear interpolation of the active `x` boundary only. The normalized semantic family converts its configured thresholds into concrete live-pass values after warmup; the absolute family uses the configured values directly. `clode.Observer.local_max` stores extrema using bounded three-sample quadratic refinement. `clode.Observer.normalized_neighborhood_return` keeps sampled anchors `x0`, then refines each stored exit time by linearly interpolating the full normalized state between the last inside sample and the first outside sample of the radius ball. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

`min_amp` is intentionally not one identical concept across all event-triggering families.

- **One-pass families** (`threshold_crossing`, `schmitt_trigger`, `local_max`): `min_amp` acts as a live-pass range gate on `event_var`.
- **Two-pass normalized threshold and Schmitt families** (`normalized_threshold_crossing`, `normalized_schmitt_trigger`): `min_amp` uses the same warmup-seeded `event_var` range that defines the normalized trigger geometry.
- **Normalized neighborhood return** (`normalized_neighborhood_return`): `min_amp` uses the warmup-seeded `event_var` range gate, while the sampled anchor threshold and normalized-ball exit geometry remain neighborhood-specific.

That difference is usually fine for steady-state workflows after a transient, but it is worth keeping in mind when you compare one-pass and two-pass event counts directly. If you need consistent amplitude semantics across a comparison, either normalize your thresholds explicitly or stick to one family.

The trigger-geometry examples now live in `examples/visualize_events_threshold_crossing.py`, `examples/visualize_events_schmitt_trigger.py`, `examples/visualize_events_localmax.py`, and `examples/visualize_events_neighborhood_return.py`. The neighborhood-return example mirrors the live sampled-anchor plus interpolated-exit geometry directly.

Recurring controls such as `min_amp`, `max_event_count`, and `max_event_timestamps` are now active across the current event-triggering family configs and remain family-local until the broader bundle-seam decision is settled.

## Observer readout inventory

The following table shows which readouts are available for each observer family:

| Readout Type | `threshold_crossing` | `normalized_threshold_crossing` | `schmitt_trigger` | `normalized_schmitt_trigger` | `local_max` | `normalized_neighborhood_return` | `summary` |
| -------------- | :----: | :----: | :----: | :----: | :----: | :----: | :----: |
| **Event timestamps** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Event values** | — | — | — | — | ✓ | — | — |
| **Event count** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Period** (time between events) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Peak count** (local maxima between events) | ✓ | ✓ | ✓ | ✓ | — | ✓ | — |
| **Schmitt duration/duty** | — | — | ✓ | ✓ | — | — | — |
| **Schmitt active dip** | — | — | ✓ | ✓ | — | — | — |
| **Amplitude** (extrema tracking) | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Trajectory max/min/mean per variable** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **Derivative extrema per variable** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **Auxiliary extrema** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

**Notes on specific readout types**:

- **Event timestamps**: Refined to sub-timestep precision with interpolation.
- **Event count**: Total number of events detected during integration. Always available for event-detector families; never available for `summary`.
- **Period**: Time between consecutive events. Meaning varies by family (inter-maxima interval for `local_max`, inter-exit time for `normalized_neighborhood_return`, etc.). Accuracy depends on event-detection strategy and oscillation regularity.
- **Amplitude**: Measured as `local_max_value - last_local_min_value` on each family's extrema-tracking channel. Threshold, Schmitt, and neighborhood-return use `feature_var`; `local_max` uses its one resolved extrema channel.
- **Peak count**: Number of local maxima detected between consecutive events. Only available for observers that track local extrema.
- **Schmitt-specific readouts**: Semantic Schmitt families now expose `up duration`, `down duration`, `duty`, and `active dip` in addition to transition times and event counts. `active dip` is the mean feature-channel value during the down state minus the last local minimum.
- **Trajectory summary statistics**: Observed extrema (`max`, `min`) and time-integrated mean for each state variable and auxiliary.

**Semantic families** (current):

- `threshold_crossing`, `normalized_threshold_crossing`: One event stream plus period/maxima/amplitude aggregates and trajectory summary statistics.
- `schmitt_trigger`, `normalized_schmitt_trigger`: Up/down transition times plus period/maxima/duration/duty/active-dip aggregates and trajectory summary statistics.
- `normalized_neighborhood_return`: Neighborhood-exit times, count, period/maxima/amplitude aggregates, sampled anchor/range fields, plus trajectory summary statistics.
- `summary`: Trajectory statistics (extrema and means) without event detection.

## Reading observer output

```python
observer_output = integrator.features()
print(observer_output.get_feature_names())
print(observer_output.get_var_mean("period"))
print(observer_output.get_event_data("up", type="time"))
```

The available names depend on the selected observer. Threshold-crossing families expose one `event` stream with `event time {index}` plus event count, period/maxima/amplitude aggregates, and trajectory summaries. `clode.Observer.local_max` exposes two extrema streams (`localmax time/value` and `localmin time/value`) plus event count and summary readouts. `clode.Observer.schmitt_trigger` and `clode.Observer.normalized_schmitt_trigger` expose separate `up transition` and `down transition` streams plus event count, period/maxima aggregates, Schmitt duration/duty readouts, the Schmitt-local `active dip` readout, and amplitude. `clode.Observer.normalized_neighborhood_return` exposes one `event` stream of interpolated normalized-ball exit times plus event count, period/maxima/amplitude aggregates, and per-variable `{var}0`/`range {var}` fields.

Trajectory-summary outputs use model variable names directly (for example `max v`, `min v`, `max dv/dt`). `clode.Observer.normalized_neighborhood_return` additionally emits per-variable center and scale fields as `{var}0` and `range {var}` (for example `v0`, `range v`).

## Custom observers

Built-in observers are supported and tested. Custom observer authoring is not part of the public API.
