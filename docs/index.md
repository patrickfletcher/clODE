# clODE documentation

clODE is a Python package for large-scale simulation of ordinary differential equation ensembles on OpenCL-capable CPUs and GPUs. This site is organized around the current supported workflows: install the package, define a model, choose a simulator, and inspect results.

## Start here

- [Install clODE](install.md) and verify that your OpenCL runtime is visible.
- [Work through the getting started guide](getting_started.md) for a first end-to-end ensemble run.
- [Inspect available platforms and devices](querying_opencl.md) before pinning platform or device IDs.

## Choose a workflow

- [Feature extraction](feature_extraction.md): compute periods, extrema, counts, and event data without storing full trajectories.
- [Trajectory simulation](trajectory_simulation.md): store time samples for plotting, post-processing, and inspection.
- [Stochastic simulation](Ornstein-Uhlenbeck_process.md): add Wiener-process terms and stochastic steppers.

## Define a model

- [Specify ODEs in Python or OpenCL](specifying_odes.md).
- [Load and convert XPP models](xpp_files.md).
- [Initialize and configure the runtime](init_runtime.md).

## Examples and reference

- [Examples](examples.md) collects runnable scripts from the repository by workflow.
- [Performance notes](performance_notes.md) documents the current benchmark scripts and the context needed to compare results responsibly.
- [Observers](observers.md) describes the built-in observer modes and related concepts.
- [API reference](api_reference.md) documents the public Python API.
- [Logging and diagnostics](logging_levels.md) covers runtime logging and PyOpenCL diagnostics.

The documentation stays focused on the current Python package and public workflows.
