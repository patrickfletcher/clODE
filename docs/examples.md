# Examples

The scripts in the `examples/` directory are a curated appendix to the workflow guides. Examples that generate plots require `matplotlib`; from a source checkout, `pip install .[docs]` is the simplest way to run the full set.

## Runtime and device inspection

- [dump_opencl_info.py](https://github.com/patrickfletcher/clODE/blob/main/examples/dump_opencl_info.py): print the OpenCL platforms and devices visible to clODE.
- [dump_device_performance.py](https://github.com/patrickfletcher/clODE/blob/main/examples/dump_device_performance.py): compare transient-solver throughput across visible devices with a Lorenz benchmark.

## Core simulation workflows

- [continuation.py](https://github.com/patrickfletcher/clODE/blob/main/examples/continuation.py): compare one long run against repeated split-window calls for both `TrajectorySimulator` and `FeatureSimulator` using the landed default continuation semantics.
- [find_steady_states.py](https://github.com/patrickfletcher/clODE/blob/main/examples/find_steady_states.py): run repeated transient solves until a small ensemble approaches steady state.
- [Ornstein_Uhlenbeck.py](https://github.com/patrickfletcher/clODE/blob/main/examples/Ornstein_Uhlenbeck.py): simulate a stochastic ensemble with the stochastic Euler stepper and compare sample statistics against the expected distribution.

## Feature extraction and event observers

- [van_der_pol_periods.py](https://github.com/patrickfletcher/clODE/blob/main/examples/van_der_pol_periods.py): measure oscillation period across an ensemble and generate the same plot used in the feature-extraction docs.
- [observe_sine_curve.py](https://github.com/patrickfletcher/clODE/blob/main/examples/observe_sine_curve.py): use a threshold observer on a simple analytic signal and inspect event timestamps.
- [spike_counting.py](https://github.com/patrickfletcher/clODE/blob/main/examples/spike_counting.py): run feature extraction across a two-parameter grid and visualize spike-count outputs.
- [visualize_events_threshold_crossing.py](https://github.com/patrickfletcher/clODE/blob/main/examples/visualize_events_threshold_crossing.py): visualize the live threshold-crossing rule shared by `threshold_crossing` and `normalized_threshold_crossing` after the normalized family derives its concrete threshold.
- [visualize_events_schmitt_trigger.py](https://github.com/patrickfletcher/clODE/blob/main/examples/visualize_events_schmitt_trigger.py): visualize the x-only semantic Schmitt state machine used by `schmitt_trigger` and `normalized_schmitt_trigger`.
- [visualize_events_localmax.py](https://github.com/patrickfletcher/clODE/blob/main/examples/visualize_events_localmax.py): visualize the sampled slope sign-change trigger and three-sample quadratic refinement used by `local_max`.
- [visualize_events_neighborhood_return.py](https://github.com/patrickfletcher/clODE/blob/main/examples/visualize_events_neighborhood_return.py): visualize the normalized neighborhood-return sampled-anchor plus interpolated normalized-ball exit geometry.
- [visualize_events_threshold2.py](https://github.com/patrickfletcher/clODE/blob/main/examples/visualize_events_threshold2.py): visualize Schmitt-trigger transitions on a two-threshold waveform.
- [visualize_events_nhood2.py](https://github.com/patrickfletcher/clODE/blob/main/examples/visualize_events_nhood2.py): visualize normalized neighborhood-return event timing in state space.

## Numerical accuracy and performance evidence

- [single_precision_accuracy.py](https://github.com/patrickfletcher/clODE/blob/main/examples/single_precision_accuracy.py): reproduce the float32 mean, time, threshold-timestamp, and local-maximum tradeoffs discussed in the numerical-accuracy docs with an inspectable NumPy mirror of the live kernel formulas.

For the underlying APIs and concepts, see [getting_started.md](getting_started.md), [feature_extraction.md](feature_extraction.md), [trajectory_simulation.md](trajectory_simulation.md), and [specifying_odes.md](specifying_odes.md). For float32 accuracy guidance, see [numerical_accuracy.md](numerical_accuracy.md). For reproducibility guidance around timings and device comparisons, see [performance_notes.md](performance_notes.md).
