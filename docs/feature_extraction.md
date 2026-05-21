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

The most directly useful current safeguards for oscillation-oriented observers are:

- `min_amp` to suppress event measurement below a chosen amplitude floor
- separate `x_up_threshold` and `x_down_threshold` values in `threshold_2` to add Schmitt-trigger-style hysteresis and reduce chatter
- `dx_up_threshold` and `dx_down_threshold` in `threshold_2` when noisy shallow crossings need an additional slope gate

Not every field is active in every built-in observer, so treat observer parameters as mode-specific rather than assuming every knob has the same effect everywhere.

Observer state also stays on the device for the duration of the solve, so its footprint depends on the observer mode, the model size, and `max_event_timestamps`. If you only need counts or summary statistics, keep `max_event_timestamps` as small as practical and prefer the lightest observer that answers the question.

`threshold_2` stores up/down transition times with inverse-linear interpolation of the active threshold boundary. When a `dx` threshold is zero, that slope gate is ignored; when it is nonzero, the stored transition time is the later of the active `x` and `dx` boundary crossings within the step. `local_max` stores extrema using bounded three-sample quadratic refinement. See [numerical_accuracy.md](numerical_accuracy.md) for empirical tradeoffs and comparisons with alternative interpolation choices.

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

## Continuation and repeated calls

`features()` continues device state by default and also continues the observer state unless
it is explicitly reinitialized. The requested `t_span` does not advance automatically.

For exact absolute-time continuation, advance the next requested window from
`get_final_time()` before calling `features()` again. This is especially important for
fixed-step methods, time-based feature accumulators, event timestamps, and non-autonomous
systems.

For autonomous systems, prefer feature windows whose local `t_span` starts near `0` when
absolute time is not part of the model. Large absolute times still coarsen stored float32
absolute timestamps even though the live observers keep elapsed-time statistics separate from
those large absolute values.

See `continuation.md` for the full continuation model and `examples/continuation.py` for a
runnable comparison between one long feature run and split-window continuation.

For empirical float32 demonstrations of amplitude floors, threshold hysteresis, derivative thresholds, compensated means, elapsed-time origins, and time accumulation limits, see [numerical_accuracy.md](numerical_accuracy.md).

## Custom observers

The built-in observers are stable and supported. Custom observer authoring is not part of the public API.
