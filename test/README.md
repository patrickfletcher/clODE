# Test Suite Layout

The supported local and CI bundles are defined in `tools/run_test_bundle.py`.
The files are not being physically moved into a final long-term directory layout yet, because the package structure is still in transition. Instead, the current suite is organized by domain through named bundles and pytest markers.

## Current domain taxonomy

- `smoke`: driver-independent packaging and frontend smoke checks.
- `frontend`: user-facing conversion and frontend integration tests.
- `runtime_api`: public runtime and simulator contract tests.
- `numerics`: exact-solution and fixed-expectation numerical regressions, including `test/core_numerics/`.
- `opencl_internal`: focused tests for internal OpenCL support layers.

## Release gating

- `smoke` runs on Linux, macOS, and Windows.
- `release` is the OpenCL-backed release gate and currently combines `frontend`, `runtime_api`, and `numerics`.
- `opencl` is kept as a compatibility alias for `release`.
- `extended` is a convenience alias for the internal OpenCL checks used for manual verification.

## Running by bundle

- `python tools/run_test_bundle.py smoke`
- `python tools/run_test_bundle.py frontend`
- `python tools/run_test_bundle.py runtime_api`
- `python tools/run_test_bundle.py numerics`
- `python tools/run_test_bundle.py release`
- `python tools/run_test_bundle.py opencl_internal`
- `python tools/run_test_bundle.py extended`

## Running by marker

- `pytest -m numerics`
- `pytest -m runtime_api`
- `pytest -m release_gate`

## Retired items

- `test/test_solver.py` was removed because it only contained skipped placeholder tests.
- `test/test_observers.py` was removed because it was a debug-only skipped script rather than an asserted regression suite.
- `test/test_clODE_utilities.py` was removed because it was a note file, not an executable test module.
- `test/test_trajectory.py` was removed because it was empty.
- `test/test_vdp_long.py` was removed because it was not part of the useful release or regression surface.
