# CI And Legacy Cleanup Plan

## Scope for the current cleanup PR

- replace the Bazel-era Linux, macOS, and Windows workflows with one clear CI workflow
- make the default CI responsibilities explicit: cross-platform smoke checks, Linux OpenCL runtime checks, docs build, and release artifact validation
- move test and docs dependencies into `pyproject.toml` extras instead of a separate `requirements.txt`
- delete non-executable placeholder tests so the remaining suite reflects real regression coverage
- keep the legacy C++ backend comparison-only and give it a smaller helper surface for local use

## Recommended CI shape

1. Cross-platform smoke jobs on Linux, macOS, and Windows should build the pure-Python package, install the built wheel, and run only driver-independent smoke tests.
2. The authoritative runtime gate should stay Linux-only and run against an explicit OpenCL ICD.
3. Docs should build in a dedicated job with `mkdocs build --strict`.
4. Release publishing should remain tag-driven and build artifacts once before publishing them.

## Optional C++ backend build options

### Option 1: Keep Bazel, add a thin helper layer

- Pros: lowest risk, no disruption to the default package, works with the current source tree.
- Cons: Bazel remains in the repository until the comparison backend is retired or rebuilt.

### Option 2: Replace Bazel with a local Makefile only

- Pros: smaller developer surface on Linux and macOS.
- Cons: poor Windows story, duplicated dependency discovery, and likely more maintenance than it saves.

### Option 3: Replace Bazel with CMake plus `scikit-build-core`

- Pros: standard cross-platform native build story, clearer Python packaging integration.
- Cons: meaningful migration work for a backend that is already comparison-only.

## Recommended path

- Do Option 1 now: keep Bazel off the default package and default CI, but hide the legacy copy/install mechanics behind a helper script and a local `clode/cpp/Makefile`.
- Only consider Option 3 if the legacy comparison backend is expected to remain supported beyond the next cleanup milestone.
- Do not spend time on Option 2 unless the comparison path becomes Linux/macOS-only by policy.

## Legacy items that can be removed now

- Bazel-specific CI workflow files
- `requirements.txt` as a duplicate dependency source
- placeholder test modules that do not execute assertions

## Legacy items that should wait until the C++ backend is retired or rebuilt

- top-level `BUILD`, `WORKSPACE`, and the `bazel/` tree
- `clode/cpp/BUILD`
- Bazel-based instructions for the standalone C++ library
- any deeper folder-structure cleanup that would break the current comparison path
