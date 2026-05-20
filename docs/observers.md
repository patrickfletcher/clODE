# Observers

Observers are the mechanism behind `FeatureSimulator`. They maintain per-ensemble state on the device and reduce trajectory information as integration proceeds, which lets clODE compute summary quantities without storing the full trajectory.

## Built-in observer modes

The current public observer modes are:

- `clode.Observer.basic`
- `clode.Observer.basic_all_variables`
- `clode.Observer.local_max`
- `clode.Observer.neighbourhood_1`
- `clode.Observer.neighbourhood_2`
- `clode.Observer.threshold_2`

Choose an observer when constructing a `FeatureSimulator`:

```python
import clode


integrator = clode.FeatureSimulator(
    src_file="test/van_der_pol_oscillator.cl",
    variables={"x": 0.0, "y": 1.0},
    parameters={"mu": 1.0},
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
)
```

## Observer parameters

Tune observer behavior with `set_observer_parameters(...)`.

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
- separated `x_up_threshold` and `x_down_threshold` values in `threshold_2` to add Schmitt-trigger-style hysteresis and reduce chatter near a boundary
- `dx_up_threshold` and `dx_down_threshold` in `threshold_2` when noisy shallow crossings need an additional slope gate

See [numerical_accuracy.md](numerical_accuracy.md) for empirical examples and practical tuning guidance.

`threshold_2` stores up/down transition times with inverse-linear interpolation of the active threshold boundary. When a `dx` threshold is zero, that slope gate is ignored; when it is nonzero, the stored transition time is the later of the active `x` and `dx` boundary crossings within the step. `local_max` stores extrema using bounded three-sample quadratic refinement. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

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
