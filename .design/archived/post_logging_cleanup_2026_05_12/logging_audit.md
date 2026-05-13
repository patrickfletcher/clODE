# Logging And Diagnostics Audit

Status: archived closeout note from the logging cleanup.
Current logging state is summarized in `.design/package_state.md` and user-facing guidance now lives in `docs/logging_levels.md`.

## Summary

- The original `clode.runtime.logging` module was compatibility state from the runtime split, not a Python-native logging design.
- The clean path is to drop the enum-style API entirely and use standard Python logging with the `clode` logger namespace.
- Explicit report helpers such as `print_opencl()` and `print_devices()` should remain explicit stdout reports, not hidden logger side effects.
- PyOpenCL already provides the diagnostics that matter most for this package: compiler-output surfacing, on-disk build caching, kernel repro capture, device characterization helpers, and pytest device parametrization.

## Current Usage Audit

| Surface | Current behavior before cleanup | Finding |
| --- | --- | --- |
| `clode/runtime/logging.py` | Stored a module-global enum and a pattern string | Compatibility shim, not a real logging backend |
| `clode.print_opencl()` | Temporarily mutates global log level and may suppress output entirely | Unpythonic for an explicit `print_*` helper |
| `OpenCLResource.print_devices()` | Returns early when `LogLevel.off` | Same issue as `print_opencl()` |
| `Simulator.print_status()` / executor status | Direct `print(...)` output | Explicit reporting is already separate from log level here |
| `LogLevel` / `set_log_level()` / `set_log_pattern()` | Exported publicly, but only carried migration-era behavior | Remove rather than preserve |
| `_opencl` runtime/build code | Emits exceptions but no standard log records | Good candidate for targeted `logging` integration |

## PyOpenCL Facilities Relevant To Diagnostics

- `PYOPENCL_COMPILER_OUTPUT`: lets PyOpenCL show compiler messages during `Program.build()` without clODE inventing a second compiler-log switch.
- `Program.build(..., cache_dir=...)`: PyOpenCL already maintains an on-disk compiler cache and allows an explicit cache directory.
- `Kernel.capture_call(...)`: writes a self-contained repro for a failing kernel launch. This is far more useful than expanding ad hoc string logging around launches.
- `pyopencl.characterize.has_src_build_cache(...)`: useful for diagnostics around whether source-build caching is expected on a given device/runtime.
- `pyopencl.characterize.get_fast_inaccurate_build_options(...)` and related helpers: useful when exposing tuning diagnostics.
- `pyopencl.tools.pytest_generate_tests_for_pyopencl(...)`: useful if clODE wants broader internal device-matrix testing later.

## Recommendation

1. Remove `LogLevel`, `get_log_level()`, `set_log_level()`, and `set_log_pattern()` from the public API.
2. Expose a minimal standard-logging surface: `configure_logging(...)` and `get_logger(...)`.
3. Treat explicit reporting helpers as explicit reporting helpers. `print_opencl()` and `simulator.print_devices()` should always print when called.
4. Route Python warnings, including `pyopencl.CompilerWarning`, through logging when requested.
5. Use targeted logger records in the internal OpenCL runtime and program-cache path instead of expanding a custom clODE logging layer.
6. Prefer PyOpenCL's own diagnostics for compiler output, cache behavior, and kernel repro capture before adding new clODE-specific knobs.

## Landed In This Pass

- `clode.runtime.logging` now exposes only `configure_logging(...)` and `get_logger(...)`.
- `configure_logging(...)` can route warnings through logging and toggle `PYOPENCL_COMPILER_OUTPUT` for subsequent builds.
- `print_opencl()` and `print_devices()` now always emit their explicit device report.
- `clode._opencl.runtime` now logs runtime selection through standard logging.
- `clode._opencl.runtime` now logs whether the selected device reports PyOpenCL source-build cache support.
- `clode._opencl.program_cache` now logs cache hits and source-build events through standard logging and includes a PyOpenCL compiler-output hint on build failure.
- Public docs now steer users toward standard logging configuration and PyOpenCL's own diagnostics knobs.

## Follow-On Recommendations

- If OpenCL build failures need more structured reproduction, add an opt-in path that uses `Kernel.capture_call(...)` rather than more string logging.
- If build-cache debugging becomes frequent, thread an optional `cache_dir` through clODE's runtime creation instead of extending the logger API.
