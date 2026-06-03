# clODE Documentation

Welcome to clODE, a Python package for large-scale ODE ensemble simulation on OpenCL-capable devices. This site is organized around actual research workflows: define a model, run an ensemble across parameter grids, and extract or inspect the results you need without always storing full trajectories. Start with the workflow pages first; the deeper single-precision and numerical-accuracy material is kept separate for when you need it.

Use this site to:

1. **Set up clODE** ([Install](install.md), [runtime setup](init_runtime.md))
2. **Learn the basics** ([Getting started](getting_started.md), [model definition](specifying_odes.md))
3. **Choose your workflow** ([Feature extraction](feature_extraction.md), [trajectory storage](trajectory_simulation.md), [stochastic simulation](Ornstein-Uhlenbeck_process.md))
4. **Find examples and reference** ([Examples](examples.md), [API reference](api_reference.md))

## Quick Navigation

**I want to...**

- Measure oscillation periods across a parameter grid → Start with [Feature extraction](feature_extraction.md)
- Store and plot time series for visualization → Start with [Trajectory simulation](trajectory_simulation.md)
- Add stochasticity to my model → See [Stochastic simulation](Ornstein-Uhlenbeck_process.md)
- Load a model from XPP/XPPAUT → See [XPP files](xpp_files.md)
- Choose the right observer for my events → See [Observers](observers.md) and the decision table in [Feature extraction](feature_extraction.md)
- Continue a solve across multiple windows → See [Continuation and repeated runs](continuation.md)
- Understand performance and single-precision trade-offs → See [Performance notes](performance_notes.md) and [Numerical accuracy](numerical_accuracy.md)
- Query my OpenCL device and platform → See [Querying OpenCL](querying_opencl.md)

## Core Concepts

- **Ensemble**: many independent ODE solves run in parallel on the device, one per work item.
- **Observer**: stateful feature detector that runs on the device during integration, accumulating periods, extrema, or event data without storing full trajectories.
- **Continuation**: repeated solves with the default `update_x0=True` continue both state and the nominal time window; use the continuation guide when you want an explicit reset or branch.

For technical reference and API details, see [API reference](api_reference.md), [Logging](logging_levels.md), and the full [Performance notes](performance_notes.md).
