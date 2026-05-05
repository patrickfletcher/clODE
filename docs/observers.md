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

## Reading observer output

```python
observer_output = integrator.features()
print(observer_output.get_feature_names())
print(observer_output.get_var_mean("period"))
print(observer_output.get_event_data("up", type="time"))
```

The available names depend on the selected observer.

## Current extension status

Built-in observers are supported and tested. User-authored custom observers are not yet a polished public API; today they still depend on internal OpenCL observer kernels and Python-side metadata definitions.
