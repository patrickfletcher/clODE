# PyOpenCL Packaging And Release Audit

## Scope

This audit checks whether the current Python packaging and release path matches the actual migration goal:

- PyOpenCL as the maintained backend
- no required Bazel or C++ build for the default Python install path
- a release process that matches semantic versioning and the real distribution shape

It also records the packaging state verified on this workspace on 2026-05-05.

## Verified Current State

- `pip wheel . --no-deps -w /tmp/clode-wheel-audit` succeeds after fixing invalid `pyproject.toml` table ordering
- `python -m build` succeeds and produces both `clode-0.9.0-cp311-cp311-linux_x86_64.whl` and `clode-0.9.0.tar.gz`
- installing the built wheel into an isolated target directory works, and `import clode` plus `clode.query_opencl()` both succeed there
- the built wheel metadata reports version `0.9.0`
- the repository tag state does not match that release identity: `git describe --tags --always --dirty` currently reports `v0.8.1-32-ged3de84-dirty`, and the newest tag is still `v0.8.1`

## Findings

### 1. The immediate packaging failure was metadata, not Bazel

The package was briefly not buildable because `pyproject.toml` placed normal `project` fields after a nested table. That made `pip wheel` fail during metadata preparation before the C++ build step even started.

That issue is now fixed locally. The current package can again build as a Bazel-backed binary wheel.

### 2. The current distribution is still a C++ and Bazel distribution

Even though the PyOpenCL backend is now implemented internally, the shipped wheel is still a platform-specific binary wheel containing:

- the compiled `clode_cpp_wrapper` extension
- the Python package
- the OpenCL kernel source files
- much of the C++ source tree and Bazel-side packaging context

This means the current package is still fundamentally distributing the legacy backend as a required part of the install story.

### 3. The public Python import path still hard-requires the C++ wrapper

The public package is not yet PyOpenCL-owned at import time.

Verified examples:

- `clode/__init__.py` imports `ProblemInfo`, `SolverParams`, and `ObserverParams` directly from `clode.cpp.clode_cpp_wrapper`
- `clode/runtime.py` imports `OpenCLResource`, runtime enums, logger types, and OpenCL query helpers from the same extension
- `clode/solver.py`, `clode/features.py`, and `clode/trajectory.py` all import wrapper structs directly
- the new `_pyopencl` implementation also still depends on C++-owned `ProblemInfo`, `SolverParams`, and `ObserverParams`

Consequence:

- switching the default backend to PyOpenCL does not remove the compiled extension from the default package
- a Bazel-free package is not possible until these public and internal type dependencies are replaced or isolated behind adapters

### 4. Versioning is currently split across incompatible sources of truth

The repository currently mixes three version identities:

- `clode.__version__ = "0.9.0"`
- `pyproject.toml` reads version dynamically from `clode.__version__`
- `setup.py` still requests `use_scm_version` through `setuptools_scm`

At the same time, git tag history still stops at `v0.8.1`.

Consequence:

- the built artifact version is not currently explained by tag history
- release automation cannot cleanly infer what should be published
- the codebase is carrying versioning machinery that is no longer acting as the authoritative source

### 5. Release automation is still aligned to the old binary-build world

The current GitHub workflows remain Bazel-centric and publish from per-OS CI jobs on push.

That is a mismatch with the desired end state.

For a true PyOpenCL-first package:

- the default Python distribution should become pure Python plus packaged `.cl` assets
- the wheel should become `py3-none-any`
- Linux wheel publication stops being a special binary-build problem because no compiled extension remains in the wheel
- publishing should be driven by release tags or an explicit release workflow, not by general push traffic from multiple platform jobs

### 6. The package currently over-ships source content and emits package-discovery warnings

The built wheel includes much more than the runtime actually needs, including the full `clode/cpp` source tree and Bazel-related files inside the package.

The sdist is broader still and currently carries most of the repository, including:

- `.github/workflows`
- `matlab/`
- `paper/`
- `test/`
- `tmp/`

The `python -m build` output also emits setuptools warnings about importable-but-undiscovered packages under:

- `clode.cpp.OpenCL`
- `clode.cpp.logging`
- `clode.cpp.observers`
- `clode.cpp.steppers`

This is acceptable while the package is still source-heavy and Bazel-backed, but it is not the right long-term layout for a clean PyOpenCL release.

### 7. Build-system hygiene still needs cleanup

The current build path still emits avoidable warnings:

- `pkg_resources` deprecation from `setup.py`
- unused or underconfigured `setuptools_scm` noise during build isolation
- deprecated license-classifier warnings from setuptools

These do not currently stop the build, but they are part of the packaging debt that should be removed before a migration-finish release.

## Recommendation

### Versioning strategy

Adopt one source of truth for the package version.

Recommended approach:

1. use tag-derived semantic versions via `setuptools_scm`
2. generate or expose the installed package version from a dedicated version module
3. remove the hardcoded `clode.__version__` literal and the current dual-source ambiguity

Recommended release semantics:

- use `0.9.0` for the current transition release line only if you need a release before the default switch
- use `0.10.0` for the first release where PyOpenCL is the default backend but a legacy C++ fallback still exists
- reserve `1.0.0` for the first release whose default install path is Bazel-free, whose runtime-support policy is explicit, and whose public API/install story you are willing to call stable

### Recommended work plan

#### PR 14: Python-Owned Public Model And Runtime Facade

Goal:

- remove the C++ wrapper from the public Python import path without changing user-facing simulator behavior

Deliverables:

- Python-owned replacements for `ProblemInfo`, `SolverParams`, and `ObserverParams`
- Python-owned runtime enums and info models where practical
- a Python `OpenCLResource`-like facade or compatibility layer for public runtime operations
- adapter code in the C++ backend that converts Python-owned models into wrapper-owned types at the leaf boundary
- `_pyopencl` internals updated to consume the Python-owned models instead of wrapper structs

Success criterion:

- `import clode` works without the C++ extension being importable, provided the PyOpenCL dependency path is installed and the C++ backend is not selected

#### PR 15: Bazel-Free Packaging And Release Transition

Goal:

- make the default Python distribution pure Python and remove Bazel from the default build path

Deliverables:

- package only the Python code and required `.cl` assets
- move kernel assets to an explicitly packaged data location rather than relying on the legacy `clode/cpp` layout
- stop shipping the compiled extension and C++ source tree in the default wheel
- quarantine the legacy C++ backend behind an explicit compatibility extra, separate package, or clearly isolated maintenance path
- replace the current push-publish workflows with a tag-driven release workflow that builds artifacts once and publishes once
- clean up the version source, `pkg_resources`, and stale `setuptools_scm` configuration debt

Success criterion:

- `python -m build` emits a `py3-none-any` wheel for the default package

#### PR 16: Default Backend Switch

Goal:

- make PyOpenCL the default backend only after the public package and release path are already PyOpenCL-owned

Deliverables:

- factory default switches to PyOpenCL
- legacy backend remains an explicit compatibility escape hatch during transition
- docs, runtime diagnostics, and support policy are updated for the default behavior

Success criterion:

- the 74-test extended reference bundle passes with the default backend on the supported runtime set

#### PR 17: Legacy C++ Extraction Or Removal

Goal:

- remove the remaining maintenance burden from the main package once transition confidence is high enough

Candidate outcomes:

- remove the C++ backend entirely
- or split it into a separate compatibility package if it still has users worth supporting

## Exit Criteria For The Migration

The migration should not be considered complete until all of the following are true:

1. `clode` imports and runs on the PyOpenCL path without requiring the C++ extension to be built or installed
2. the default package builds without Bazel and publishes a pure-Python wheel
3. the default backend is PyOpenCL on the supported runtime set
4. the release version is derived from one authoritative source and matches git tags
5. the release workflow publishes from a dedicated release/tag path rather than from general CI pushes
