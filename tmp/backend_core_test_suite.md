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

Today this combined gate is 20 tests and is the suite that should stay green while the backend seam lands.

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

This file is intentionally small. Its job is to lock down the current backend contract, not to become a second broad API suite.

## Runtime Guidance

- Keep the authoritative gate single precision for now.
- Use the public simulator classes, but keep the test intent backend-neutral.
- Prefer exact assertions where the reference formulas make that realistic.
- Keep tolerances file-local and explicit where adaptive or stochastic behavior requires them.

On this Linux workspace, the stable local command is:

```bash
CLODE_TEST_PLATFORM_ID=1 CLODE_TEST_DEVICE_ID=0 /home/fletcherpa/envs/clode/bin/python -m pytest test/core_numerics test/test_backend_contracts.py -q
```

The environment-variable override lives in `test/core_numerics/helpers.py` so the tests do not hardcode local device IDs.

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
