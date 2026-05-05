# Minimal Backend-Phase Test Suite

## Purpose

This document defines the authoritative fast test surface to protect clODE while the backend seam and PyOpenCL migration work proceed.

It has two layers:

- a numerical-core suite that treats the Python layer as a thin harness for launching kernels and reading buffers
- a small backend-contract suite that locks down wrapper and runtime behaviors the next phase depends on

Everything else under `test/` is reference material, historical coverage, or out of scope for the initial backend migration gate.

## Authoritative Test Files

### Numerical core authority

- `test/core_numerics/test_transient.py`
- `test/core_numerics/test_trajectory.py`
- `test/core_numerics/test_features_basicall.py`
- `test/core_numerics/test_stochastic.py`

### Backend contract authority

- `test/test_backend_contracts.py`
- `test/test_backend_rhs_source.py`
- `test/test_pyopencl_models.py`
- `test/test_pyopencl_source_builder.py`
- `test/test_pyopencl_runtime.py`
- `test/test_pyopencl_buffers.py`
- `test/test_pyopencl_structs.py`
- `test/test_pyopencl_transient_backend.py`
- `test/test_pyopencl_trajectory_backend.py`
- `test/test_pyopencl_feature_backend.py`

Today this combined gate is 54 tests and is the suite that should stay green through the subsequent backend-migration phases.

Current implementation state:

- the public simulators now construct their execution backend through `clode/_backends/factory.py`
- the active reference path is `clode/_backends/cpp.py`, which wraps the current pybind C++ runtime
- public simulators now prepare an internal `RhsSource` object with source text and digest before backend construction
- the internal `clode/_pyopencl/` package now contains the immutable build, source, and program models plus backend-specific error types
- the internal `clode/_pyopencl/` package now also contains a static kernel registry and deterministic phase-one source builder
- the internal `clode/_pyopencl/` package now also contains a real PyOpenCL runtime selector and runtime-scoped program cache
- the internal `clode/_pyopencl/` package now also contains a common-state buffer manager that centralizes the current flatten and reshape rules
- the internal `clode/_pyopencl/` package now also contains device-matched host struct helpers for `SolverParams` and `ObserverParams`
- the internal `clode/_pyopencl/` package now also contains a transient executor that can be selected internally through `_CLODE_BACKEND=pyopencl`
- the internal `clode/_pyopencl/` package now also contains trajectory buffer support and a trajectory executor behind the same internal backend selector
- the internal `clode/_pyopencl/` package now also contains observer metadata and layout sizing, feature buffers, and a feature executor behind the same internal backend selector

## Canonical Models

The numerical suite is built around four committed OpenCL fixtures.

| Fixture | Why it exists | Main coverage |
| --- | --- | --- |
| `stable_linear.cl` | exact deterministic reference with closed-form state and derivative | transient exactness, continuation, fixed-step timing |
| `stable_linear_aux.cl` | exact deterministic reference with auxiliary outputs | trajectory aux coverage, `basicall` aux summaries |
| `ornstein_uhlenbeck.cl` | minimal stochastic reference with known stationary moments | seeded reproducibility and ensemble statistics |
| `hopf_normal_form.cl` | exact on-cycle oscillatory reference | adaptive trajectory and `basicall` cycle summaries |

No test in the authoritative suite should rely on Python-generated RHS source when a committed `.cl` fixture can express the behavior under test.

## Numerical Core Coverage

### `test_transient.py`

This file currently covers:

- RK4 stable-linear final-state exactness
- RK4 deterministic ensemble exactness
- Dormand-Prince stable-linear final-state exactness
- deterministic split-window continuation parity for RK4 and Dormand-Prince

The transient gate is intentionally about final states, final times, and continuation-sensitive state evolution. It does not try to exhaustively test wrapper setters.

### `test_trajectory.py`

This file currently covers:

- RK4 stable-linear `t`, `x`, and `dx` samples against exact references
- RK4 stable-linear-aux exact auxiliary samples
- fixed-step storage contract for `nout` and `max_store`
- Dormand-Prince Hopf trajectory parity against the exact on-cycle solution at returned sample times

### `test_features_basicall.py`

This file currently covers:

- RK4 `basicall` state summaries on a settled linear window
- RK4 `basicall` auxiliary summaries on a settled linear-aux window
- Dormand-Prince `basicall` Hopf cycle summaries on the limit cycle
- fixed-step step-count behavior through the observer path

The current gate deliberately uses settled windows for the decaying linear cases rather than startup-transient windows.

### `test_stochastic.py`

This file currently covers:

- seeded stochastic Euler reproducibility on the same backend
- long-run Ornstein-Uhlenbeck mean and variance against theory

The stochastic gate does not currently include a continuation test. That was intentionally kept out of the minimal gate to preserve reliability and runtime budget.

## Backend Contract Coverage

`test/test_backend_contracts.py` protects the non-numerical behaviors that matter before backend work starts:

- `get_tspan()` returns the current device window instead of `None`
- cached `get_final_state()` calls are safe on repeat access
- `set_solver_parameters()` resets the live device `dt` buffer on an existing simulator
- observer switching rebuilds and runs through the public `FeatureSimulator` path
- `features(initialize_observer=...)` refreshes cached feature results instead of returning stale data
- `max_event_timestamps` changes rebuild the compiled feature program with the expected `N_STORE_EVENTS` value

`test/test_backend_rhs_source.py` protects the PR 4 source-preparation path:

- file-backed RHS inputs are prepared as source text plus digest
- Python-callable RHS inputs are prepared as source text plus digest
- source digest changes when source text changes

`test/test_pyopencl_models.py` protects the PR 5 model and diagnostics layer:

- `ProblemShape` mirrors solver dimensions from `ProblemInfo`
- `BuildKey` remains hashable and enforces the observer and stored-event contract
- `ProgramBundle` enforces build-key and kernel-handle consistency
- PyOpenCL error types retain source, options, and other diagnostic context

`test/test_pyopencl_source_builder.py` protects the PR 6 source-preparation layer:

- the Python registry matches the current C++ stepper and observer define map
- transient, trajectory, and feature source assembly use the expected entrypoint order
- build options preserve the current compile-time specialization knobs
- the build key changes when RHS source changes while kernel-tree digest stays stable
- conservative RHS validation rejects source that does not define `getRHS(...)`

`test/test_pyopencl_runtime.py` protects the PR 7 runtime and compile layer:

- explicit runtime selection chooses one requested OpenCL device
- transient source bundles compile successfully through PyOpenCL on the configured test device
- the program cache is runtime-scoped and reuses cache hits by `BuildKey`
- build failures preserve source text, build options, and build log details

`test/test_pyopencl_buffers.py` protects the PR 8 buffer layer:

- `ArrayLayout` preserves the current Fortran-order flatten and reshape contract
- common transient buffers allocate the expected byte sizes for state, parameters, RNG state, `dt`, and `tf`
- upload and download helpers round-trip problem data, `dt`, and RNG state through PyOpenCL buffers
- solver-parameter packing matches the current C-style struct layout used by the kernels

`test/test_pyopencl_structs.py` protects the struct-handling rules established after the PyOpenCL audit:

- `SolverParams` and `ObserverParams` use device-matched PyOpenCL dtypes rather than relying on handwritten host assumptions
- packing helpers preserve field values across single and double precision paths
- the simplest double-precision observer-state layout demonstrates why raw byte-count formulas are unsafe for later feature-backend work

`test/test_pyopencl_transient_backend.py` protects the PR 9 transient execution path:

- deterministic transient parity is checked against the current C++ backend
- continuation of `x0`, `dt`, and final time matches the current backend contract
- seeded stochastic ensembles are reproducible on the PyOpenCL backend and remain numerically aligned with the C++ path
- updating solver parameters still rewrites the live `dt` buffer on the PyOpenCL path

`test/test_pyopencl_trajectory_backend.py` protects the PR 10 trajectory execution path:

- deterministic trajectory samples match the current C++ backend for `t`, `x`, and `dx`
- auxiliary trajectory samples remain aligned with the current C++ backend
- `nout`, `max_store`, and `n_stored` semantics are pinned on the PyOpenCL path
- trajectory storage reallocates correctly when `max_store` changes without changing the public `TrajectoryOutput` behavior

`test/test_pyopencl_feature_backend.py` protects the PR 11 feature execution path:

- `basicall` feature statistics remain numerically aligned with the current C++ backend
- the `localmax` rebuild contract updates `N_STORE_EVENTS` and feature-name exposure correctly on the PyOpenCL path
- a two-pass observer path (`thresh2`) remains aligned with the current C++ backend

Feature backend edge cases that remain active after PR 11:

- `ObserverData` is observer-specific opaque state and must be sized from an explicit layout model rather than copied from the current C++ host byte-count formulas
- feature execution uses a distinct `initializeObserver` kernel before the `features` kernel when observer state needs to be initialized or reinitialized
- two-pass observers such as `nhood2` and `thresh2` rely on that initialize pass to compute detector thresholds before feature collection begins
- the zero-parameter Python-callable RHS constructor edge case is still a current C++-path limitation and stays out of scope for the PyOpenCL parity work
- the Intel CPU OpenCL runtime on this workspace remains unstable for the `localmax` feature rebuild path, so milestone validation currently prefers the NVIDIA runtime

This file is intentionally small. Its job is to lock down the current backend contract, not to become a second broad API suite.

## Runtime Guidance

- Keep the authoritative gate single precision for now.
- Use the public simulator classes, but keep the test intent backend-neutral.
- Prefer exact assertions where the reference formulas make that realistic.
- Keep tolerances file-local and explicit where adaptive or stochastic behavior requires them.

On this Linux workspace, the stable local command is:

```bash
CLODE_TEST_PLATFORM_ID=0 CLODE_TEST_DEVICE_ID=0 /home/fletcherpa/envs/clode/bin/python -m pytest test/core_numerics/test_transient.py test/core_numerics/test_trajectory.py test/core_numerics/test_features_basicall.py test/core_numerics/test_stochastic.py test/test_backend_contracts.py test/test_backend_rhs_source.py test/test_pyopencl_models.py test/test_pyopencl_source_builder.py test/test_pyopencl_runtime.py test/test_pyopencl_buffers.py test/test_pyopencl_structs.py test/test_pyopencl_transient_backend.py test/test_pyopencl_trajectory_backend.py test/test_pyopencl_feature_backend.py -q
```

For the broader PR 11 acceptance bundle on this workspace, add `test/test_vdp.py`, `test/test_features.py`, `test/test_aux_values.py`, and `test/test_ornl_thompson_a1.py` to the same command. That expanded bundle is currently 63 tests and passed on the stable NVIDIA runtime.

The environment-variable override lives in `test/core_numerics/helpers.py` so the tests do not hardcode local device IDs.

Important runtime-selection note:

- `clinfo -l` reports Intel as platform `0` and NVIDIA as platform `1` on this machine
- the stable backend-validation command above currently targets the NVIDIA runtime with `CLODE_TEST_PLATFORM_ID=0` and `CLODE_TEST_DEVICE_ID=0`
- do not assume that `clinfo -l` platform ordering matches the runtime-selection order seen by the current backend validation path on this workspace

## Out Of Scope For This Gate

The following remain outside the authoritative backend-phase gate:

- `test/test_function_converter.py`
- `test/test_xpp_parser.py`
- `test/test_runtime.py`
- `test/test_logger.py`
- `test/test_opencl_builtins.py`
- `test/test_features.py`
- `test/test_aux_values.py`
- `test/test_vdp.py`
- `test/test_ornl_thompson_a1.py`
- `test/test_solver.py`
- `test/test_trajectory.py`
- `test/test_observers.py`

Those files may still be useful as reference or later expansion points, but they are not the pinned authority for the migration work.
