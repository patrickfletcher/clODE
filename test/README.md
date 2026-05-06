# Test Suite Layout

The supported CI bundles are defined in `tools/run_test_bundle.py`.

## Cross-platform smoke bundle

- Runs on Linux, macOS, and Windows.
- Covers package-safe, driver-independent checks.
- Includes `test/test_function_converter.py`.
- Includes `test/test_pyopencl_models.py`.
- Includes `test/test_pyopencl_source_builder.py`.

## Linux OpenCL bundle

- Runs on Linux with an installed OpenCL ICD.
- This is the default runtime gate for pull requests.
- Includes `test/core_numerics/`.
- Includes `test/test_backend_contracts.py`.
- Includes `test/test_backend_rhs_source.py`.
- Includes `test/test_pyopencl_runtime.py`.
- Includes `test/test_runtime.py`.
- Includes `test/test_logger.py`.

## Extended bundle

- Kept for manual or scheduled runtime verification.
- Includes broader runtime/API regressions such as `test/test_vdp.py`, `test/test_xpp_parser.py`, and the `test/test_pyopencl_*backend.py` modules.
- Not part of the default pull-request gate.

## Long bundle

- Contains opt-in long-running checks.
- Currently includes `test/test_vdp_long.py`.

## Retired legacy placeholders

- `test/test_solver.py` was removed because it only contained skipped placeholder tests.
- `test/test_observers.py` was removed because it was a debug-only skipped script rather than an asserted regression suite.
- `test/test_clODE_utilities.py` was removed because it was a note file, not an executable test module.
