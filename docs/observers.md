# Observers

Observers are the mechanism behind `FeatureSimulator`. They maintain per-ensemble state on the device and reduce trajectory information as integration proceeds, which lets clODE compute summary quantities without storing the full trajectory.

## Built-in observer modes

The current public observer modes are:

- `clode.Observer.summary`
- `clode.Observer.basic`
- `clode.Observer.basic_all_variables`
- `clode.Observer.local_max`
- `clode.Observer.threshold_crossing`
- `clode.Observer.normalized_threshold_crossing`
- `clode.Observer.schmitt_trigger`
- `clode.Observer.normalized_schmitt_trigger`
- `clode.Observer.threshold_2`
- `clode.Observer.neighborhood_return`
- `clode.Observer.neighbourhood_1`
- `clode.Observer.neighbourhood_2`

`clode.Observer.local_max` is the canonical extrema observer (with `clode.Observer.local_extremum` retained as a compatibility alias). `clode.Observer.threshold_2` and `clode.Observer.neighbourhood_2` remain retained legacy full-featured observers when you need their broader readout bundles.

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
- `clode.Observer.threshold_2`: retained legacy fully featured normalized Schmitt trigger when you need the broader period, duty, active-dip, and trajectory-summary readout bundle

Threshold crossing and Schmitt triggering stay separate user-facing concepts even though equal up/down thresholds can collapse to a degenerate single-boundary case internally. They expose different readout bundles and answer different workflow questions.

## Choosing an extremum or neighborhood workflow

- `clode.Observer.local_max`: canonical maxima-oriented workflow with IMI, amplitude, all-state summaries, and separate local-maximum/local-minimum event streams.
- `clode.Observer.neighborhood_return`: lean two-pass normalized neighborhood-return detector. Warmup fixes the normalization range, the live pass anchors on the first sampled point where `event_var` drops below `anchor_threshold`, and the observer stores linearly refined neighborhood-exit times plus event count.
- `clode.Observer.neighbourhood_2`: retained legacy fully featured normalized neighborhood-return workflow when you also need period, peak, and broader summary outputs. It now shares the same sampled-anchor plus interpolated normalized-ball exit timing as the lean semantic family.

## Observer parameters

Threshold-family observers now have preferred semantic config objects. Use `observer_configuration=` at construction time or `set_observer_configuration(...)` when you want a family-specific knob set instead of the broad compatibility bundle.

Use `clode.ThresholdCrossingConfig` for both threshold-crossing families:

```python
integrator = clode.FeatureSimulator(
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

`x_up_threshold` must be greater than or equal to `x_down_threshold`. Equality is allowed and collapses the hysteresis band to one shared boundary while keeping the Schmitt-style up/down readout streams. If you need derivative-gated crossings, stay on the legacy `clode.Observer.threshold_2` compatibility surface.

Use `clode.LocalMaximumConfig` with `clode.Observer.local_max` for maxima-triggered extrema tracking:

```python
integrator.set_observer_configuration(
    clode.LocalMaximumConfig(
        event_var="x",
        max_event_timestamps=16,
    ),
    observer=clode.Observer.local_max,
)
```

Use `clode.NeighborhoodReturnConfig` with `clode.Observer.neighborhood_return` when you want the lean normalized neighborhood-return trigger:

```python
integrator.set_observer_configuration(
    clode.NeighborhoodReturnConfig(
        event_var="x",
        anchor_threshold=0.25,
        radius=0.15,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.neighborhood_return,
)
```

`anchor_threshold` picks the falling `event_var` section that latches the anchor point `x0`, while `radius` sets the size of the normalized full-state neighborhood around that anchor. This can be a good fit when a one-variable threshold crossing is too ambiguous to identify one recurrence cleanly. On a simple limit cycle it behaves like a light-weight return-map trigger around one anchored point on the orbit; on a more complex or multi-lobed cycle it can still distinguish nearby passes because the exit test uses the full normalized state rather than only the scalar `event_var`.

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

Not every observer field affects every built-in observer. For the current threshold families, the most directly useful safeguards are:

- `min_amp` for suppressing oscillation/event measurements below a chosen amplitude floor
- `threshold` plus `direction` in `clode.ThresholdCrossingConfig` for the threshold-crossing families, with the selected observer deciding whether that threshold is absolute or warmup-derived
- separate `x_up_threshold` and `x_down_threshold` values in `clode.SchmittTriggerConfig` for the Schmitt families, with the selected observer deciding whether those values are absolute or warmup-derived fractions
- the retained legacy `clode.Observer.threshold_2` compatibility surface when noisy shallow crossings need derivative gates in addition to the value thresholds

See [numerical_accuracy.md](numerical_accuracy.md) for empirical examples and practical tuning guidance.

`clode.Observer.schmitt_trigger` and `clode.Observer.normalized_schmitt_trigger` store up/down transition times with inverse-linear interpolation of the active `x` boundary only. The normalized semantic family converts its configured thresholds into concrete live-pass values after warmup; the absolute family uses the configured values directly. `clode.Observer.threshold_2` keeps the heavier legacy period, duty, and active-state readout bundle and still supports derivative gates on the compatibility surface. `clode.Observer.local_max` stores extrema using bounded three-sample quadratic refinement. `clode.Observer.neighborhood_return` and `clode.Observer.neighbourhood_2` keep sampled anchors `x0`, but now refine each stored exit time by linearly interpolating the full normalized state between the last inside sample and the first outside sample of the radius ball. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

`min_amp` is intentionally not one identical concept across all threshold-style families. 

- **One-pass absolute families** (`threshold_crossing`, `schmitt_trigger`): `min_amp` acts as a live-pass range gate on `event_var`. The threshold is compared directly against `min_amp` to suppress events from oscillations below that amplitude.
- **Two-pass normalized families** (`normalized_threshold_crossing`, `normalized_schmitt_trigger`): `min_amp` is compared against the warmup-derived amplitude used to define the normalized thresholds. This means the effective amplitude floor adapts to your warmup trajectory.
- **Legacy `threshold_2`**: Follows the two-pass semantics (warmup-derived amplitude floor).

That difference is usually fine for steady-state workflows after a transient, but it is worth keeping in mind when you compare one-pass and two-pass event counts directly. If you need consistent amplitude semantics across a comparison, either normalize your thresholds explicitly or stick to one family.

The trigger-geometry examples now live in `examples/visualize_events_threshold_crossing.py`, `examples/visualize_events_schmitt_trigger.py`, `examples/visualize_events_localmax.py`, and `examples/visualize_events_neighborhood_return.py`. The neighborhood-return example mirrors the live sampled-anchor plus interpolated-exit geometry directly. The older `threshold_2`, `local_max`, and `neighbourhood_2` visualizations remain useful legacy comparisons.

Recurring controls such as `min_amp`, `max_event_count`, and `max_event_timestamps` still stay on the family configs where they are active rather than moving into a separate shared oscillation bundle. The threshold-crossing readout-selection surface is a narrower pilot toward the longer-term goal of letting users choose specific readout subsets per observer family.

## Observer readout inventory

The following table shows which readouts are available for each observer family:

| Readout Type | `threshold_crossing` | `normalized_threshold_crossing` | `schmitt_trigger` | `normalized_schmitt_trigger` | `local_max` | `neighborhood_return` | `neighbourhood_2` (legacy) | `summary` |
|--------------|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|:----:|
| **Event timestamps** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Event values** | — | — | — | — | ✓ | — | — | — |
| **Event count** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Period** (time between events) | — | — | — | — | ✓ | ✓ | ✓ | — |
| **Amplitude** (extrema tracking) | — | — | — | — | ✓ | — | — | — |
| **Peak count** (maxima per period) | — | — | — | — | — | ✓ | ✓ | — |
| **Up duration** (Schmitt-specific) | — | — | — | — | — | — | — | — |
| **Down duration** (Schmitt-specific) | — | — | — | — | — | — | — | — |
| **Duty cycle** (Schmitt-specific) | — | — | — | — | — | — | — | — |
| **Trajectory max/min/mean per variable** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |
| **Derivative extrema per variable** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | — |
| **Auxiliary extrema** | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ | ✓ |

**Notes on specific readout types**:

- **Event timestamps**: Refined to sub-timestep precision with interpolation.
- **Event count**: Total number of events detected during integration. Always available for event-detector families; never available for `summary`.
- **Period**: Time between consecutive events. Meaning varies by family (inter-maxima interval for `local_max`, inter-exit time for `neighbourhood_2`, etc.). Documented per family below. Accuracy depends on event-detection strategy and oscillation regularity.
- **Amplitude**: Measured as `local_max_value - last_local_min_value` in the variable specified by `feature_var`. Only available for extrema-tracking families.
- **Peak count**: Number of local maxima detected between consecutive events. Only available for observers that track local extrema.
- **Schmitt-specific readouts** (up/down duration, duty cycle): Only available for Schmitt trigger families. Currently only `threshold_2` (legacy) computes these; semantic Schmitt families store transition times but not the per-period statistics.
- **Trajectory summary statistics**: Observed extrema (`max`, `min`) and time-integrated mean for each state variable and auxiliary. Not yet available for semantic extremum families; current only on legacy and summary observers. Future work: add to all event-detector families for consistency.

**Semantic families** (lean, preferred):
- `threshold_crossing`, `normalized_threshold_crossing`: Narrow event-stream readouts (timestamps ± count). Minimal state footprint.
- `schmitt_trigger`, `normalized_schmitt_trigger`: Up/down transition times ± count. Minimal state footprint.
- `neighborhood_return`: Neighborhood-exit times, count. Minimal state footprint.
- `summary`: Trajectory statistics (extrema and means) without event detection.

**Legacy families** (retained, full-featured):
- `threshold_2`: Normalized Schmitt with period, duty, peak counts, and trajectory summaries. Heavier state footprint.
- `local_max`: Local-extremum detection with period, amplitude, peak counts, and trajectory summaries. Heavier state footprint.
- `neighbourhood_2`: Normalized neighborhood-return with period, peak counts, and trajectory summaries. Heavier state footprint.

## Reading observer output

```python
observer_output = integrator.features()
print(observer_output.get_feature_names())
print(observer_output.get_var_mean("period"))
print(observer_output.get_event_data("up", type="time"))
```

The available names depend on the selected observer. Threshold-crossing families expose one `event` stream with `event time {index}` plus count. `clode.Observer.local_max` exposes two extrema streams (`local maximum time/value` and `local minimum time/value`) plus event count and summary readouts. `clode.Observer.schmitt_trigger` and `clode.Observer.normalized_schmitt_trigger` expose separate `up transition` and `down transition` streams plus an event count, while `clode.Observer.threshold_2` keeps those streams and also retains the heavier legacy period, duty, and active-state readouts. `clode.Observer.neighborhood_return` exposes one `event` stream of interpolated normalized-ball exit times, while `clode.Observer.neighbourhood_2` retains the broader legacy periodicity bundle around the same sampled-anchor plus interpolated-exit trigger.

## Custom observers

Built-in observers are supported and tested. Custom observer authoring is not part of the public API.
