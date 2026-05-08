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

- `clode/__init__.py`: curated export surface and package version.
- `clode/problem/*`: problem definition, source loading, Python/XPP ingestion, and OpenCL equation authoring helpers.
- `clode/runtime/*`: public OpenCL query, device selection, and logging helpers.
- `clode/simulation/*`: solver params, simulator orchestration, and result containers.
- `clode/observers/*`: observer enums, parameter schema, and feature-name catalog helpers.
- Root compatibility barrels remain for public flat imports such as `solver.py`, `features.py`, `trajectory.py`, `function_converter.py`, `xpp_parser.py`, and `opencl_builtins.py`.

### Internal OpenCL implementation

- `clode/_opencl/models.py`: build keys, source bundles, problem shape, and precision enums.
- `clode/_opencl/runtime.py`: explicit single-device context and queue creation.
- `clode/_opencl/registry.py`: stepper and observer define registry plus entrypoint mapping.
- `clode/_opencl/source_builder.py`: kernel assembly, build options, and kernel-tree digesting.
- `clode/_opencl/program_cache.py`: runtime-scoped OpenCL program cache.
- `clode/_opencl/structs.py`: device-matched struct dtypes.
- `clode/_opencl/buffers.py`: buffer allocation plus flatten/reshape rules.
- `clode/_opencl/observer_metadata.py`: runtime-specific observer-data struct modeling.
- `clode/_opencl/executors.py`: transient, trajectory, and feature executors.

The historical backend-shim and binding-named compatibility layers have been removed.

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
- `test/test_simulation_contracts.py`: current simulation and observer behavior contracts.
- `test/test_problem_rhs_source.py`: RHS source-ingestion and digest coverage.
- `test/test_opencl_models.py`, `test/test_opencl_source_builder.py`, `test/test_opencl_runtime.py`, `test/test_opencl_buffers.py`, `test/test_opencl_structs.py`: canonical internal OpenCL support-layer tests.
- `tools/probe_opencl_runtime.py`: distinguishes runtime/compiler failures from clODE kernel failures.
- `docs/init_runtime.md`: current runtime-selection story, including the still-public `device_ids` note.

## Packaging Quick Lookup

- `pyproject.toml`: authoritative packaging and dependency configuration.
- `MANIFEST.in`: current sdist include and prune rules.
- Kernel assets ship as package data.
- Docs and tests are not required for runtime execution. Tests are still useful for downstream verification, while docs are the easier thing to omit if sdist slimming becomes desirable.
