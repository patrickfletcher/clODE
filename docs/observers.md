# Observers

Observers are the mechanism behind `FeatureSimulator`. They maintain per-ensemble state on the device and reduce trajectory information as integration proceeds, which lets clODE compute summary quantities without storing the full trajectory.

## Built-in observer modes

The current public observer modes are:

- `clode.Observer.summary`
- `clode.Observer.basic`
- `clode.Observer.basic_all_variables`
- `clode.Observer.local_max`
- `clode.Observer.threshold_crossing`
- `clode.Observer.neighbourhood_1`
- `clode.Observer.neighbourhood_2`
- `clode.Observer.schmitt_trigger`

The older pass-count names `clode.Observer.threshold_1` and `clode.Observer.threshold_2` remain supported as compatibility aliases.

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

## Observer parameters

Threshold-family observers now have preferred semantic config objects. Use `observer_configuration=` at construction time or `set_observer_configuration(...)` when you want a family-specific knob set instead of the broad compatibility bundle.

```python
integrator = clode.FeatureSimulator(
    ...,
    observer_configuration=clode.ThresholdCrossingConfig(
        event_var="x",
        threshold=-40.0,
        direction=clode.EventDirection.rising,
        max_event_timestamps=16,
    ),
)
```

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
    )
)
```

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

Not every observer field affects every built-in observer. For the current built-in set, the most directly useful numerical safeguards are:

- `min_amp` for suppressing oscillation/event measurements below a chosen amplitude floor
- separated `x_up_threshold` and `x_down_threshold` values in `clode.Observer.schmitt_trigger` to add Schmitt-trigger-style hysteresis and reduce chatter near a boundary
- `dx_up_threshold` and `dx_down_threshold` in `clode.Observer.schmitt_trigger` when noisy shallow crossings need an additional slope gate

See [numerical_accuracy.md](numerical_accuracy.md) for empirical examples and practical tuning guidance.

`clode.Observer.schmitt_trigger` stores up/down transition times with inverse-linear interpolation of the active threshold boundary. When a `dx` threshold is zero, that slope gate is ignored; when it is nonzero, the stored transition time is the later of the active `x` and `dx` boundary crossings within the step. `local_max` stores extrema using bounded three-sample quadratic refinement. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

For threshold-trigger workflows, the current observer catalog is easier to reason about if you keep two ideas separate:

- `clode.Observer.threshold_crossing` is a one-pass directional threshold crossing in absolute units of `event_var`. Use it when the event variable is expected to cross a physically meaningful value such as a voltage level.
- `clode.ThresholdCrossingConfig` is the preferred way to configure that family because it exposes `threshold` and `direction` directly instead of the broader compatibility names.
- `clode.Observer.schmitt_trigger` is a two-pass warmup-derived fractional Schmitt trigger. Its `x_up_threshold` and `x_down_threshold` inputs are interpreted as fractions of the observed warmup amplitude of `event_var`, so it is often more robust across parameter sweeps where the event-variable range changes.
- `clode.SchmittTriggerConfig` is the preferred way to configure that family because it exposes only the currently relevant Schmitt-trigger knobs, including `min_amp` and the optional `dx` gates.
- `dx_up_threshold` and `dx_down_threshold` in `clode.Observer.schmitt_trigger` are optional extra gates that are most helpful on noisy or stochastic traces.
- A Schmitt trigger with equal up/down thresholds can degenerate to a simple threshold crossing internally, but treating threshold crossing and Schmitt triggering as separate public concepts keeps configuration and output expectations clearer.

## Reading observer output

```python
observer_output = integrator.features()
print(observer_output.get_feature_names())
print(observer_output.get_var_mean("period"))
print(observer_output.get_event_data("up", type="time"))
```

The available names depend on the selected observer.

## Custom observers

Built-in observers are supported and tested. Custom observer authoring is not part of the public API.
