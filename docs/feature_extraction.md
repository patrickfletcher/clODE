# Feature Extraction

`FeatureSimulator` computes trajectory statistics on the device while integrating, which keeps memory use low even for large ensembles. Instead of storing every time sample, it runs a stateful **observer** that accumulates features such as periods, extrema, event counts, or event timestamps, and returns only those summaries when the solve completes.

This is the right tool when your research question is about ensemble behavior or aggregate properties rather than individual trajectories. You gain memory efficiency, speed, and the ability to scale to very large ensembles.

## Which observer should I use?

clODE includes seven observer families, each optimized for a different workflow. Pick the one that matches your feature question:

| Your Question | Best Observer | Configuration | Output |
| --- | --- | --- | --- |
| What are the min/max/mean of state variables or their time derivatives? | `Observer.summary` | `SummaryObserverSelection` (optional) | Summary statistics, no event detection |
| When do I cross a fixed threshold (e.g., voltage = 0)? | `Observer.threshold_crossing` | `ThresholdCrossingConfig` with absolute threshold | Event times, event count, period/peak-count/amplitude aggregates |
| When do I cross a threshold _relative to my amplitude range_? | `Observer.normalized_threshold_crossing` | `ThresholdCrossingConfig` with fractional threshold | Event times (warmup-scaled), event count, period/peak-count/amplitude aggregates |
| When do I transition across two state boundaries (hysteresis)? | `Observer.schmitt_trigger` | `SchmittTriggerConfig` with up/down thresholds | Up/down transition times, event count, period/peak-count/duration/duty/amplitude aggregates |
| When do I transition with hysteresis _relative to amplitude_? | `Observer.normalized_schmitt_trigger` | `SchmittTriggerConfig` (warmup-scaled) | Up/down transition times (warmup-scaled), event count, period/peak-count/duration/duty/amplitude aggregates |
| When do I reach local maxima or minima? | `Observer.local_max` | `LocalMaximumConfig` | Max/min times and values, IMI, amplitude, trajectory summaries |
| When do I return to a sampled neighborhood on a limit cycle? | `Observer.normalized_neighborhood_return` | `NeighborhoodReturnConfig` | Exit times, event count, period/peak-count/amplitude aggregates, sampled anchor state |

Use `get_feature_names()` on a simulator or the returned `ObserverOutput` to see all available feature names for your chosen observer.

## Configuring Your Chosen Observer

Each observer family has a preferred semantic configuration class. Pass it to the constructor or call `set_observer_configuration(...)` before running `features()`:

### Threshold-Crossing Observers

Use `ThresholdCrossingConfig` for either `threshold_crossing` family:

```python
simulator = clode.FeatureSimulator(
    ...,
    observer=clode.Observer.normalized_threshold_crossing,
    observer_configuration=clode.ThresholdCrossingConfig(
        event_var="x",
        feature_var="x",
        threshold=0.75,  # fraction of warmup amplitude for normalized_threshold_crossing
        direction=clode.EventDirection.rising,
        max_event_timestamps=16,
    ),
)
```

- `Observer.threshold_crossing`: `threshold` is an absolute value in state-variable units.
- `Observer.normalized_threshold_crossing`: `threshold` is a fraction of the warmup-pass amplitude of `event_var`.

`event_var` chooses the threshold geometry and `min_amp` gate for the threshold families. `feature_var` chooses the extrema/amplitude channel, so threshold and Schmitt now share the same trigger-versus-measurement split when you want event detection on one state variable and oscillation readouts on another.

### Schmitt-Trigger Observers

Use `SchmittTriggerConfig` for either `schmitt_trigger` family:

```python
simulator.set_observer_configuration(
    clode.SchmittTriggerConfig(
        event_var="x",
        x_up_threshold=0.3,
        x_down_threshold=0.2,  # must be ≤ x_up_threshold
        max_event_timestamps=16,
    ),
    observer=clode.Observer.normalized_schmitt_trigger,
)
```

- `Observer.schmitt_trigger`: thresholds are absolute values.
- `Observer.normalized_schmitt_trigger`: thresholds are warmup-derived fractions of the `event_var` amplitude.

### Local Maximum Observer

Use `LocalMaximumConfig`:

```python
simulator.set_observer_configuration(
    clode.LocalMaximumConfig(
        event_var="x",
        max_event_timestamps=16,
    ),
    observer=clode.Observer.local_max,
)
```

### Neighborhood Return Observer

Use `NeighborhoodReturnConfig`:

```python
simulator.set_observer_configuration(
    clode.NeighborhoodReturnConfig(
        event_var="x",
        feature_var="x",
        anchor_threshold=0.25,
        radius=0.15,
        max_event_timestamps=16,
    ),
    observer=clode.Observer.normalized_neighborhood_return,
)
```

`anchor_threshold` still picks the falling `event_var` section that latches the sampled anchor. `feature_var` selects the extrema/amplitude channel for neighborhood-return readouts; keep it aligned with `event_var` for one-channel behavior, or point it at another state variable when the return trigger and the measured oscillation should differ.

### Summary Observer (Optional Config)

Use `SummaryObserverSelection` to choose which summaries to compute:

```python
selection = clode.SummaryObserverSelection(
    state={"x": ("max", "mean"), "y": ("min",)},
    slope={"x": "min"},
)

simulator = clode.FeatureSimulator(
    ...,
    observer=clode.Observer.summary,
    summary_selection=selection,
)
```

If you do not provide a summary selection, `Observer.summary` returns mean, min, and max for all state variables.

## Reading Observer Results

`features()` returns an `ObserverOutput` object with the computed feature data. Common accessors are:

```python
output = simulator.features()

# Get feature names
names = output.get_feature_names()

# Get summary statistics
mean_val = output.get_var_mean("x")
min_val = output.get_var_min("x")
max_val = output.get_var_max("x")

# Get time-derivative statistics (slope)
min_slope = output.get_var_min_slope("x")

# Get event data
event_times = output.get_event_data("event", type="time")
event_count = output.get_event_data("event", type="count")
```

Each observer family exposes different feature names and streams. For example:

- `Observer.summary` exposes state variable summaries and optionally derivative summaries.
- `Observer.threshold_crossing` exposes one `event` stream with crossing times.
- `Observer.schmitt_trigger` exposes separate `up` and `down` streams with transition times.
- `Observer.local_max` exposes `local maximum` and `local minimum` streams plus IMI and amplitude.

Use `get_feature_names()` to explore what is available for your observer configuration, or see the examples below for visual demonstrations.

## Continuation and Repeated Calls

By default, `features()` continues device state and observer state from where the previous call left off. The requested `t_span` does not auto-advance.

For multi-window feature extraction:

- Before calling `features()` again, advance the next window from `get_final_time()` instead of using the requested end time from the previous call.
- This is especially important for event timestamps and non-autonomous systems.

For full continuation semantics, solver diagnostics, and split-window examples, see [continuation.md](continuation.md).
