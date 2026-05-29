# Performance notes

clODE is designed for large ensembles of independent ODE solves and for workflows where feature extraction can stay on the device instead of storing full trajectories. This page is intentionally conservative: it documents what can be reproduced from the repository today without turning incomplete benchmark work into headline claims.

This page is about throughput. For reproducible demonstrations of single-precision numerical error accumulation and the current mitigation tradeoffs, see [numerical_accuracy.md](numerical_accuracy.md).

## Current status

- Formal cross-package benchmark claims are still in progress.
- The repository already includes scripts that let you compare clODE configurations and visible OpenCL devices on your own hardware.
- Results depend strongly on the OpenCL runtime, device, precision, stepper, observer mode, ensemble size, and whether trajectories or only final states/features are stored.

## Reproducible scripts in the repository

- `examples/dump_opencl_info.py`: record the OpenCL platforms and devices visible to clODE before running any comparison.
- `examples/dump_device_performance.py`: run a Lorenz-system transient benchmark across visible devices using single-precision RK4 and report min, median, and max times across repeated runs.
- `examples/single_precision_accuracy.py`: reproduce the float32 mean-accumulation and time-accumulation demonstrations discussed in the numerical-accuracy docs.

When publishing or sharing a result, record the exact script, any local edits, and the commit you ran.

## What to report with a benchmark

At minimum, include:

- operating system
- Python version
- clODE version or commit
- PyOpenCL version
- OpenCL platform, device, and driver/runtime details
- model and model size
- stepper and precision
- `dt`, tolerances, and `t_span`
- observer mode or trajectory-storage mode
- ensemble size and repetition count

Without this context, comparisons are usually not meaningful.

## Practical guidance

- Warm up the target device before timing repeated runs.
- Compare like with like: same model, same precision, same stepper, same output mode.
- Treat transient-only, feature-extraction, and full-trajectory workloads as different performance regimes.
- Treat observer mode and `max_event_timestamps` as part of the workload definition. Event-driven observers such as `normalized_schmitt_trigger`, `normalized_neighborhood_return`, and `local_max` maintain more persistent per-instance state than `summary`, and retained event timestamps scale linearly with `max_event_timestamps`.
- If you only need summary features, keep event timestamp retention small and prefer the lightest observer that answers the question.
- Prefer median timing across repeated runs over a single best-case measurement.
- Record whether the runtime is CPU-backed, GPU-backed, or PoCL-backed.

## What is still missing

- a stable benchmark matrix that covers transient, feature, and trajectory workloads
- cross-package comparisons presented with enough methodology to survive review
- benchmark plots and narratives suitable for the paper and package landing pages

For now, use this page as a reproducibility note and use the example scripts as the supported starting point for local performance investigations.
