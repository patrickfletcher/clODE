# Observers

Observers are the mechanism behind `FeatureSimulator`. They maintain per-ensemble state on the device and reduce trajectory information as integration proceeds, which lets clODE compute summary quantities without storing the full trajectory.

## Built-in observer modes

The current public observer modes are:

- `clode.Observer.summary`
- `clode.Observer.basic`
- `clode.Observer.basic_all_variables`
- `clode.Observer.local_extremum`
- `clode.Observer.local_max`
- `clode.Observer.threshold_crossing`
- `clode.Observer.normalized_threshold_crossing`
- `clode.Observer.schmitt_trigger`
- `clode.Observer.normalized_schmitt_trigger`
- `clode.Observer.threshold_2`
- `clode.Observer.neighborhood_return`
- `clode.Observer.neighbourhood_1`
- `clode.Observer.neighbourhood_2`

`clode.Observer.local_extremum` and `clode.Observer.neighborhood_return` are the preferred semantic observers for lean extremum and normalized neighborhood-return workflows. `clode.Observer.local_max`, `clode.Observer.threshold_2`, and `clode.Observer.neighbourhood_2` remain retained legacy full-featured observers when you need their broader readout bundles.

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
- `clode.Observer.normalized_schmitt_trigger`: two-pass Schmitt trigger where the `x_*` and optional `dx_*` thresholds are interpreted as warmup-derived fractions
- `clode.Observer.threshold_2`: retained legacy fully featured normalized Schmitt trigger when you need the broader period, duty, active-dip, and trajectory-summary readout bundle

Threshold crossing and Schmitt triggering stay separate user-facing concepts even though equal up/down thresholds can collapse to a degenerate single-boundary case internally. They expose different readout bundles and answer different workflow questions.

## Choosing an extremum or neighborhood workflow

- `clode.Observer.local_extremum`: lean local-extremum detector on one variable. `clode.LocalExtremumConfig.polarity` selects maxima, minima, or both, and the readout keeps event times, event values, and event count.
- `clode.Observer.local_max`: retained legacy maxima-oriented workflow when you still need IMI, amplitude, all-state summaries, and separate `localmax` and `localmin` event streams.
- `clode.Observer.neighborhood_return`: lean two-pass normalized neighborhood-return detector. Warmup fixes the normalization range, the live pass anchors on the first point where `event_var` drops below `anchor_threshold`, and the observer stores neighborhood-exit times plus event count.
- `clode.Observer.neighbourhood_2`: retained legacy fully featured normalized neighborhood-return workflow when you also need period, peak, and broader summary outputs.

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

Use `clode.SchmittTriggerConfig` for both semantic Schmitt families and the retained `clode.Observer.threshold_2` workflow:

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
        dx_up_threshold=0.0,
        dx_down_threshold=0.0,
    ),
    observer=clode.Observer.normalized_schmitt_trigger,
)
```

Use `clode.LocalExtremumConfig` with `clode.Observer.local_extremum` when you want a polarity-selectable lean extremum stream:

```python
integrator.set_observer_configuration(
    clode.LocalExtremumConfig(
        variable="x",
        polarity=clode.ExtremumPolarity.maximum,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.local_extremum,
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

When you call `set_observer_configuration(...)` with one of the shared config classes, pass `observer=` whenever you are switching families or when the simulator was not already constructed with the desired observer. The config class alone does not determine absolute versus warmup-derived interpretation.

Single-observer configs such as `clode.LocalExtremumConfig` and `clode.NeighborhoodReturnConfig` are unambiguous, so clODE can infer their observer when you construct a simulator from the config alone. Passing `observer=` explicitly is still useful when you want the selection to be obvious at the call site.

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
- separate `x_up_threshold` and `x_down_threshold` values in `clode.SchmittTriggerConfig` for the Schmitt families, with the selected observer deciding whether those values are absolute or warmup-derived
- `dx_up_threshold` and `dx_down_threshold` in `clode.SchmittTriggerConfig` when noisy shallow crossings need an additional slope gate; `dx_down_threshold` is specified as a positive magnitude for the required negative downward slope

See [numerical_accuracy.md](numerical_accuracy.md) for empirical examples and practical tuning guidance.

`clode.Observer.schmitt_trigger`, `clode.Observer.normalized_schmitt_trigger`, and `clode.Observer.threshold_2` store up/down transition times with inverse-linear interpolation of the active threshold boundary. The normalized semantic family and `threshold_2` convert their configured thresholds into concrete live-pass values after warmup; the absolute family uses the configured values directly. When a `dx` threshold is zero, that slope gate is ignored; when it is nonzero, the stored transition time is the later of the active `x` and `dx` boundary crossings within the step. `threshold_2` additionally keeps the heavier legacy period, duty, and active-state readout bundle. `clode.Observer.local_extremum` and `clode.Observer.local_max` store extrema using bounded three-sample quadratic refinement. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

Recurring controls such as `min_amp`, `max_event_count`, and `max_event_timestamps` currently stay on the family configs where they are active rather than moving into a separate shared oscillation bundle.

## Reading observer output

```python
observer_output = integrator.features()
print(observer_output.get_feature_names())
print(observer_output.get_var_mean("period"))
print(observer_output.get_event_data("up", type="time"))
```

The available names depend on the selected observer. The threshold-crossing families expose one `threshold` event stream. `clode.Observer.local_extremum` exposes one `local_extremum` event stream with timestamps and event values, while `clode.Observer.local_max` keeps separate `localmax` and `localmin` streams plus the heavier legacy summaries. `clode.Observer.schmitt_trigger` and `clode.Observer.normalized_schmitt_trigger` expose separate `up` and `down` streams plus an event count, while `clode.Observer.threshold_2` keeps those streams and also retains the heavier legacy period, duty, and active-state readouts. `clode.Observer.neighborhood_return` exposes one `neighborhood_return` stream, while `clode.Observer.neighbourhood_2` retains the broader legacy periodicity bundle around the same general trigger.

## Custom observers

Built-in observers are supported and tested. Custom observer authoring is not part of the public API.
