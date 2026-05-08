# clODE Package State

## Snapshot

- clODE is now a pure-Python, PyOpenCL-only package.
- The public API is centered on `Simulator`, `TrajectorySimulator`, and `FeatureSimulator`.
- Runtime-critical OpenCL source assets live under `clode/kernels/` and ship as package data.
- Historical migration and bug notes are archived under `tmp/archived/`; they remain useful for rationale, but they are not the source of truth for the live package.
- The current package still exposes a few transition-era surfaces, especially `device_ids` and the coupling between solver configuration and output/storage capacity.

## Top-Level Repo Map

| Path | Role | Notes |
| --- | --- | --- |
| `clode/` | Maintained package source | All runtime-critical code lives here |
| `docs/` | MkDocs source | Current user-facing Python docs |
| `examples/` | Current example scripts and sample RHS files | Live Python/OpenCL examples |
| `test/` | Supported regression suite | Organized by bundle and marker rather than final directory layout |
| `tools/` | Development and diagnostics helpers | Includes the test-bundle runner and OpenCL probe |
| `tmp/archived/` | Archived migration, cleanup, and bug notes | High-signal design history, but some statements are historical only |
| `dist/`, `site/`, `clode.egg-info/` | Generated packaging or docs artifacts | Useful locally, not authoritative source |
| `paper/` | Paper describing this project (target: Journal of Open Source Software) | Orthogonal to the runtime and package internals |

## Package Map

### Public layer

- `clode/__init__.py`: export surface and package version.
- `clode/runtime.py`: public OpenCL query and selection helpers; PyOpenCL-only backend resolution.
- `clode/types.py`: Python-owned `ProblemInfo`, `SolverParams`, and `ObserverParams`.
- `clode/solver.py`: base ensemble simulator; owns source preparation, ensemble shaping, cache invalidation, and transient execution flow.
- `clode/trajectory.py`: trajectory storage/output wrapper.
- `clode/features.py`: observer-backed feature extraction wrapper.
- `clode/function_converter.py` and `clode/xpp_parser.py`: RHS generation and conversion tools.
- `clode/opencl_builtins.py`: Python names that map to OpenCL builtins for equation authoring.

### Internal seam

- `clode/_backends/protocol.py`: internal simulator, trajectory, and feature backend interfaces.
- `clode/_backends/factory.py`: backend creation and runtime-selection normalization.
- `clode/_backends/rhs.py`: `RhsSource` text-plus-digest model.

This seam is now thin. It mostly isolates the public wrappers from the PyOpenCL executors, but it still reflects the transition-era architecture.

### PyOpenCL implementation

- `clode/_pyopencl/models.py`: build keys, source bundles, problem shape, and precision enums.
- `clode/_pyopencl/runtime.py`: explicit single-device context and queue creation.
- `clode/_pyopencl/registry.py`: stepper and observer define registry plus entrypoint mapping.
- `clode/_pyopencl/source_builder.py`: kernel assembly, build options, and kernel-tree digesting.
- `clode/_pyopencl/program_cache.py`: runtime-scoped OpenCL program cache.
- `clode/_pyopencl/structs.py`: device-matched struct dtypes via PyOpenCL.
- `clode/_pyopencl/buffers.py`: buffer allocation plus flatten/reshape rules.
- `clode/_pyopencl/observer_metadata.py`: observer feature names and `ObserverData` layout modeling.
- `clode/_pyopencl/executors.py`: transient, trajectory, and feature backends.

### Kernel tree

- `clode/kernels/transient.cl`: base transient kernel.
- `clode/kernels/trajectory.cl`: trajectory storage path.
- `clode/kernels/initializeObserver.cl` and `clode/kernels/features.cl`: feature and observer lifecycle.
- `clode/kernels/steppers/*.clh`: explicit and adaptive stepper implementations.
- `clode/kernels/observers/*.clh`: built-in observer implementations.
- `clode/kernels/odedriver.cl`: currently unused unified-driver prototype; intentionally kept as deferred future-design context, not active runtime code.

## Runtime Model

- One OpenCL work-item advances one ODE instance.
- Problem arrays are flattened in Fortran order on the Python side.
- Kernel builds remain compile-time specialized by precision, stepper, observer, and problem dimensions.
- `FeatureSimulator` uses persistent observer state plus an optional warmup kernel for two-pass observers.
- Continuation is object-stateful: runs normally advance `x0`, keep device `dt`, keep RNG state, and for features keep observer state unless reinitialized.
- Absolute-time continuation is still not a first-class internal model; callers manage `t_span` explicitly.

## Current Constraints And Live Debt

- The runtime is effectively single-device only, even though `device_ids` is still accepted in public constructors.
- `SolverParams` still mixes integration controls with trajectory-storage controls (`max_store`, `nout`).
- `ObserverParams` defaults are still duplicated between the public constructors and the dataclass.
- Observer metadata is explicit but still hardcoded through large conditional logic rather than a cleaner observer-definition model.
- Optional event storage still participates in compile-time observer layout and buffer sizing.
- Split-window continuation is still numerically inconsistent for at least:
  - `basicall` feature continuation
  - seeded stochastic Euler continuation
- `tmp/archived/` is useful for rationale and bug archaeology, but some archived statements about the old wrapper path are now historical only.

## Test And Tooling Quick Lookup

- `tools/run_test_bundle.py`: authoritative bundle map.
- `test/core_numerics/`: exact-solution and kernel-level regression backbone.
- `test/test_backend_contracts.py`: current behavior contracts that were unstable during migration.
- `tools/probe_opencl_runtime.py`: distinguishes runtime/compiler failures from clODE kernel failures.
- `docs/init_runtime.md`: current runtime-selection story, including the still-public `device_ids` note.

## Packaging Quick Lookup

- `pyproject.toml`: authoritative packaging and dependency configuration.
- `MANIFEST.in`: current sdist include and prune rules.
- Kernel assets ship as package data.
- Docs and tests are not required for runtime execution. Tests are still useful for downstream verification, while docs are the easier thing to omit if sdist slimming becomes desirable.
