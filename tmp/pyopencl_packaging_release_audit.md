# PyOpenCL Packaging And Release Audit

## Scope

This audit checks whether the current Python packaging and release path matches the actual migration goal:

- PyOpenCL as the maintained backend
- no required Bazel or C++ build for the default Python install path
- a release process that matches semantic versioning and the real distribution shape

It also records the packaging and release state verified on this workspace through PR 16.

## Verified Current State

- `python -m build` now succeeds through the `pyproject.toml` path without invoking Bazel and emits `clode-0.10.0-py3-none-any.whl` plus `clode-0.10.0.tar.gz`
- the default wheel contains Python modules plus packaged OpenCL assets under `clode/kernels/`; it no longer ships the wrapper extension or the `clode/cpp` source tree
- an isolated install of the built wheel works outside the repository checkout, and the installed package resolves to the PyOpenCL runtime path when the wrapper is absent
- focused PR 15 regression coverage passes on this workspace: `test/test_pyopencl_source_builder.py` and `test/test_pyopencl_runtime.py` both pass against the packaged kernel root
- `twine check` passes for the built artifacts
- `matlab/`, `samples/`, `paper/`, `tmp/`, and `.github/` are now outside the default source-distribution payload
- the current branch version is now `0.10.0`, matching the intended first pure-Python default-backend release line; git tags still need to be advanced when the release is cut

## Findings

### 1. PR 15 converted the default package to a pure-Python endpoint

The default build path now runs entirely through `pyproject.toml` and a stub `setup.py`.

Consequences:

- `python -m build` no longer compiles the wrapper extension
- the default wheel is now `py3-none-any`
- Bazel is no longer part of the default Python release path

### 2. Runtime assets are now packaged explicitly instead of piggybacking on `clode/cpp`

The OpenCL source files needed at runtime now live under `clode/kernels/` for packaging purposes.

Consequences:

- the wheel contains the actual runtime assets and not the legacy host-side source tree
- the C++ compatibility path can still consume the same relative kernel layout when pointed at the packaged root

### 3. The package now defaults to PyOpenCL even in a source checkout

In an installed wheel where the wrapper is absent, backend resolution uses PyOpenCL by default.

In a source checkout that still contains a locally built wrapper binary, backend resolution now still uses PyOpenCL by default and only selects the legacy path when `_CLODE_BACKEND=cpp` is set explicitly.

Consequence:

- the selector policy is now aligned with the packaging story rather than depending on local checkout artifacts

### 4. Release automation is now structurally aligned, but version policy is not yet settled

Push CI jobs no longer publish to PyPI. Release publication now lives in one dedicated tag-driven workflow.

Remaining mismatch:

- the package version is now sourced consistently from `clode.__version__`, but the repository tags still lag that version and no final semver policy has been adopted for the migration finish line

### 5. Package artifacts are materially leaner, though the sdist still carries docs, examples, and tests

The default source distribution now excludes several clearly non-runtime trees, including `.github/`, `paper/`, and `tmp/`, in addition to the already-pruned stale `matlab/` and `samples/` trees.

That is acceptable for the current transition stage. Further sdist slimming is optional rather than a blocker.

### 6. The modern PyOpenCL-only endpoint is now concrete rather than aspirational

For the desired end state, the package should converge on standard modern Python packaging practices rather than a customized binary-extension build path.

That target should look like this:

- `pyproject.toml` is the authoritative packaging configuration
- the default wheel is pure Python and built with `python -m build`
- package data is declared explicitly for runtime `.cl` assets
- build-time behavior does not depend on `setup.py` side effects
- release publishing is tag-driven and uses one dedicated workflow
- CI separates validation from publishing
- package version comes from one authoritative semver-compatible source

That means the migration finish line is not just "PyOpenCL works". It is also "the package looks like a normal modern Python package".

## Recommendation

### Versioning strategy

Adopt one source of truth for the package version.

Recommended approach:

1. decide whether `clode.__version__` remains the long-term source of truth or whether a later tag-derived workflow is worth reintroducing
2. make that choice explicit in the release workflow and release notes before the first post-PR15 publication
3. ensure git tags match the published artifact version before any release claiming the pure-Python packaging transition

Recommended release semantics:

- use `0.9.0` for the current transition release line only if you need a release before the default switch
- use `0.10.0` for the first release where PyOpenCL is the default backend but a legacy C++ fallback still exists
- reserve `1.0.0` for the first release whose default install path is Bazel-free, whose runtime-support policy is explicit, and whose public API/install story you are willing to call stable

### Recommended work plan

### Modern packaging target for the PyOpenCL-only endpoint

Packaging:

1. move package configuration ownership fully into `pyproject.toml`
2. remove the default compiled extension from `ext_modules`
3. package only Python modules plus the OpenCL source assets needed at runtime
4. stop shipping stale or non-runtime trees in sdists and wheels, including `matlab/` and `samples/`
5. remove `pkg_resources` and stale `setup.py` compatibility logic from the default build path
6. unify version sourcing, ideally via a dedicated version module populated from release tags

CI and release workflows:

1. split validation workflows from publishing workflows
2. stop publishing from general `push` jobs on Linux, macOS, and Windows
3. run the normal test matrix on push and pull request events
4. run build-artifact verification in one packaging workflow using `python -m build`
5. publish only from tagged releases or a dedicated release dispatch workflow
6. use trusted publishing or a single explicit PyPI publish job instead of repeating upload logic across OS jobs

Repository scope:

1. keep stale maintenance surfaces out of the release path
2. treat `matlab/` and `samples/` as out of scope for Python package artifacts
3. revisit whether they should remain in the repository at all after the backend migration settles

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

Status update:

- Landed on this workspace: public model types and runtime compatibility objects are now Python-owned, `clode/_backends/factory.py` no longer imports the C++ adapter eagerly, and a subprocess regression test now blocks `clode.cpp` imports while verifying the PyOpenCL public path still imports and constructs a simulator

#### PR 15: Bazel-Free Packaging And Release Transition

Goal:

- make the default Python distribution pure Python and remove Bazel from the default build path

Deliverables:

- package only the Python code and required `.cl` assets
- move kernel assets to an explicitly packaged data location rather than relying on the legacy `clode/cpp` layout
- stop shipping the compiled extension and C++ source tree in the default wheel
- stop shipping stale repository trees in package artifacts, including `matlab/` and `samples/`
- quarantine the legacy C++ backend behind an explicit compatibility extra, separate package, or clearly isolated maintenance path
- replace the current push-publish workflows with a tag-driven release workflow that builds artifacts once and publishes once
- clean up the version source, `pkg_resources`, and stale `setuptools_scm` configuration debt

Success criterion:

- `python -m build` emits a `py3-none-any` wheel for the default package

Status update:

- landed on this workspace
- validated by `python -m build`, `twine check`, focused PyOpenCL runtime/source-builder tests, and an isolated wheel smoke test that resolves to PyOpenCL without the wrapper

#### PR 16: Default Backend Switch

Goal:

- make PyOpenCL the default backend only after the public package and release path are already PyOpenCL-owned

Deliverables:

- factory default switches to PyOpenCL
- legacy backend remains an explicit compatibility escape hatch during transition
- docs, runtime diagnostics, and support policy are updated for the default behavior

Success criterion:

- the extended reference bundle passes with the default backend on the supported runtime set

Status update:

- landed on this workspace
- the selector now defaults to PyOpenCL even when a locally built wrapper binary exists in `clode/cpp/`
- legacy comparison remains possible through explicit `_CLODE_BACKEND=cpp` selection after building the wrapper in a source checkout

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
