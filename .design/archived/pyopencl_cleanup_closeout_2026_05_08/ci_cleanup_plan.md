# CI, Testing, And Packaging Status

## Implemented in the current cleanup pass

- Bazel is no longer part of the default Python package path or CI path.
- The remaining Bazel support is trimmed to the legacy wrapper build only.
- The old Bazel Python packaging layer is removed.
- `pyproject.toml` is now the single source of truth for package, docs, and test dependencies.
- The supported Python floor is aligned with the actual `pyopencl` dependency requirement: Python 3.10+.
- The current tests are grouped by domain through named bundles and pytest markers, without forcing a premature final directory reorganization.

## Minimal Bazel surface that remains justified

The retained Bazel pieces are only those needed to build `//clode/cpp:clode_cpp_wrapper` and the optional standalone C++ library targets under `clode/cpp/`.

### Retained

- `WORKSPACE`
- `bazel/clode_http_archive.bzl`
- `bazel/repositories.bzl`
- `bazel/external_deps.bzl`
- `bazel/repository_locations.bzl`
- `bazel/repository_locations_utils.bzl`
- `bazel/python/**`
- `bazel/remote_config/**`
- `bazel/external/fmtlib.BUILD`
- `bazel/external/spdlog.BUILD`
- `bazel/external/pybind11.BUILD`
- `bazel/external/opencl_windows.BUILD`
- `clode/cpp/**/BUILD`

### Removed as dead packaging support

- top-level `BUILD`
- `clode/BUILD`
- `bazel/get_python_libs.py`
- `bazel/external/python.BUILD`
- unused `rules_python` and `rules_apple` setup from `WORKSPACE`
- unused external Python tarball metadata in the Bazel repository definitions

## Testing recommendation for the current transition phase

Do not do a large physical test-tree migration yet.

The better near-term move is:

1. keep the authoritative kernel-level numerical tests in `test/core_numerics/`
2. keep higher-level exact-solution and observer-output regressions in place
3. expose a clearer domain taxonomy through bundle names and pytest markers
4. reserve a future directory rewrite for after the final package structure settles

That is now the implemented approach.

## Current test taxonomy

- `smoke`: driver-independent packaging and frontend smoke checks
- `frontend`: function conversion, OpenCL builtins, and XPP conversion checks
- `runtime_api`: public runtime, backend-selection, and simulator contract checks
- `numerics`: exact-solution, observer-output, and `core_numerics` regressions
- `pyopencl_internal`: focused tests for the internal PyOpenCL support layers
- `legacy_cpp_comparison`: tests that still compare against the optional C++ wrapper backend
- `release`: the current OpenCL-backed release gate, combining `frontend`, `runtime_api`, and `numerics`

## CI audit against PyOpenCL upstream

Upstream `inducer/pyopencl` CI focuses on runtime diversity more than Python-version fan-out. In particular, it uses:

- Linux POCL as a predictable default runtime gate
- separate Intel OpenCL jobs where they still want vendor-specific coverage
- cross-platform wheel builds as a separate concern from runtime CI

That points to the following conclusions for clODE:

1. using `pocl-opencl-icd` for the default GitHub-hosted Linux runtime gate is the right choice
2. Intel OpenCL is not a good default PR gate on GitHub runners because setup is more brittle and the signal is less predictable
3. several Python versions do make sense for clODE, but mostly for packaging and import smoke, not for the expensive runtime job

## Recommended CI shape going forward

- Cross-platform smoke: Linux, macOS, and Windows on a small Python matrix that includes the supported minimum and a current upper version.
- Default runtime gate: Linux plus POCL on one current Python version.
- Docs: strict MkDocs build on pull requests.
- Release: tag-driven artifact build and publish, with build happening once.

## Packaging audit

For the current pure-Python package, `setuptools` remains a reasonable backend.

The current best-practice adjustments are:

- use `pyproject.toml` as the primary configuration file
- avoid a `setup.py` shim when it is no longer needed
- avoid duplicate dependency declarations in a separate `requirements.txt`
- align `requires-python` with the real dependency floor
- avoid unexplained upper bounds on the build backend unless there is a known breakage to pin around

Those adjustments are now reflected in the repository.

## Recommended next step after this pass

Once the final package layout is settled and the legacy C++ wrapper is either removed or formally frozen, do a second-stage test-tree rewrite that moves files physically into domain-oriented directories. Until then, the current bundle-plus-marker taxonomy is the lower-risk path.
