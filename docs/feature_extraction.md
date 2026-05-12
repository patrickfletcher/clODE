# Feature extraction

`FeatureSimulator` computes trajectory statistics on the device while integrating, which keeps memory use low even for large ensembles. This is the right tool when you want summary quantities such as periods, extrema, event counts, or event timestamps instead of full trajectories.

## Built-in observers

clODE currently ships the following observer modes through the public `Observer` enum:

- `Observer.basic`: basic summary statistics for one variable
- `Observer.basic_all_variables`: basic summary statistics for all state variables
- `Observer.local_max`: local-maximum event tracking
- `Observer.neighbourhood_1`
- `Observer.neighbourhood_2`
- `Observer.threshold_2`: threshold-based event tracking and period-style measurements

The exact feature names depend on the observer. Use `get_feature_names()` on a configured simulator or `get_feature_names()` on the resulting `ObserverOutput` to inspect what is available.

## Example

The example below measures the period of the Van der Pol oscillator across an ensemble of `mu` values.

This example is also stored as [examples/van_der_pol_periods.py](https://github.com/patrickfletcher/clODE/blob/main/examples/van_der_pol_periods.py) so the docs and the runnable script stay aligned.

```py source run
--8<-- "examples/van_der_pol_periods.py"
```

## Configuring an observer

Observer configuration lives in `ObserverParams` and can be updated through `set_observer_parameters(...)`.

Common options include:

- `event_var`: which variable is used for event detection
- `feature_var`: which variable is used for feature readout
- `max_event_count`: how many events to accumulate
- `max_event_timestamps`: how many event timestamps to retain
- `min_amp`, `min_imi`, `nhood_radius`, `x_up_threshold`, `x_down_threshold`, `dx_up_threshold`, `dx_down_threshold`, `eps_dx`

Example:

```python
integrator.set_observer_parameters(
    event_var="x",
    feature_var="x",
    max_event_timestamps=16,
    min_amp=0.1,
)
```

## Reading results

`features()` returns an `ObserverOutput` object. Common accessors are:

- `get_feature_names()`
- `get_var_mean(name)`
- `get_var_min(name)`
- `get_var_max(name)`
- `get_var_count(name)`
- `get_event_data(name, type="time")`

## Current customization status

The built-in observers are stable and supported. Authoring completely custom observers is still an internal workflow tied to the OpenCL observer kernels and Python-side observer metadata. That is an active design area for the post-migration cleanup, but it is not yet a polished public extension API.
