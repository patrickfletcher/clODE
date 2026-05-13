# Backend Overhaul PR Task Breakdown

## Purpose

This document turns the migration design into PR-sized work items with checkpoints and acceptance criteria.

It is intentionally more granular than the phase plan in `.design/pyopencl_backend_design.md`. The design phases are the roadmap. This document is the execution plan.

## Current Progress

- Completed: PR 0 scope and test lock
- Completed: PR 1 numerical-core and backend-contract test lock
- Completed: PR 2 backend protocol and factory
- Completed: PR 3 C++ backend adapter and public-wrapper rewiring
- Completed: PR 4 RHS source object integration
- Completed: PR 5 PyOpenCL core models and errors
- Completed: PR 6 PyOpenCL registry and source builder
- Completed: PR 7 PyOpenCL runtime and program cache
- Completed: PR 8 Buffer manager for common state
- Completed: PR 9 PyOpenCL transient backend
- Completed: PR 10 trajectory backend
- Completed: PR 11 feature backend
- Completed: PR 12 extended reference suite and rollout guardrails
- Completed: PR 13 dependency surfacing and runtime diagnostics
- Completed: PR 14 Python-owned public types and runtime facade
- Completed: PR 15 Bazel-free packaging and release transition
- Completed: PR 16 default backend switch
- Current milestone audit: the authoritative backend gate, the extended PyOpenCL reference bundle, the default-backend `test/test_vdp.py` smoke, and the transition example at `examples/pyopencl_ornstein_uhlenbeck.py` are the current validation targets on the stable NVIDIA runtime on this workspace; packaging and release-state notes are recorded in `.design/pyopencl_packaging_release_audit.md`
- Struct-handling audit result: host-populated OpenCL structs now use device-matched PyOpenCL dtypes for `SolverParams` and `ObserverParams`; `ObserverData` remains a deferred feature-backend concern and must not reuse the legacy byte-count formulas
- Trajectory note: the current kernel contract counts `max_store` as total storage slots including the initial sample at slot 0; this existing behavior was respected in the new regression coverage and was not changed here
- Deferred follow-up: a zero-parameter Python-callable RHS can still trip a current C++ backend construction-time edge case on this Linux workspace; keep it documented but out of scope for the current PR sequence
- Feature note for PR 11: feature execution is a two-kernel lifecycle (`initializeObserver` plus `features`) and observer continuation depends on preserving the opaque per-ensemble `ObserverData` buffer across calls
- Runtime note: `clinfo -l` reports Intel as platform `0` and NVIDIA as platform `1` on this workspace, but the current stable backend-validation command uses `CLODE_TEST_PLATFORM_ID=0` and `CLODE_TEST_DEVICE_ID=0` to target the NVIDIA runtime; the Intel CPU runtime remains unstable for the `localmax` rebuild path

## Planning Rules

- Each PR should have one primary deliverable.
- Each PR should preserve existing public behavior unless the acceptance criteria explicitly say otherwise.
- Each PR must update docs when it changes authoritative scope, tests, or rollout state.
- No PR should mix backend seam work, PyOpenCL runtime work, and post-parity cleanup.

## Milestone Overview

| Milestone | Intent | Status |
| --- | --- | --- |
| M0 | Lock scope and tests | Complete |
| M1 | Land backend seam | Complete |
| M2 | Land PyOpenCL build primitives | Complete |
| M3 | Land transient parity | Complete |
| M4 | Land trajectory parity | Complete |
| M5 | Land feature parity | Complete |
| M6 | Switch defaults and simplify later | Not started |

## PR 0: Scope And Test Lock

Priority: P0

Goal:

- freeze the authoritative scope and pinned test suite before implementation starts

Deliverables:

- `.design/archived/backend_migration_history_2026_05_07/backend_overhaul_scope_map.md`
- `.design/archived/backend_migration_history_2026_05_07/backend_core_test_suite.md`
- `.design/backend_pr_task_breakdown.md`

Checkpoint:

- maintainers agree which tests are authoritative and which are excluded

Acceptance criteria:

- all three docs exist and are internally consistent
- excluded placeholder files are named explicitly
- the pinned suite is narrower than the full `test/` tree

## PR 1: Add Missing Contract Tests

Priority: P0

Goal:

- fill the current coverage gaps before backend changes start

Deliverables:

- `test/test_backend_contracts.py`
- `test/core_numerics/`

Checkpoint:

- the new tests run only against the current C++ backend initially

Acceptance criteria:

- new tests are added exactly as scoped in `.design/archived/backend_migration_history_2026_05_07/backend_core_test_suite.md`
- no backend seam changes yet
- the updated core suite passes on the current C++ backend

## PR 2: Backend Protocol And Factory

Priority: P0

Status:

- Complete

Goal:

- introduce the internal backend boundary without changing behavior

Deliverables:

- `clode/_backends/protocol.py`
- `clode/_backends/factory.py`
- minimal `clode/_backends/__init__.py`

Checkpoint:

- public modules can import the protocol and factory without using them yet everywhere

Acceptance criteria:

- the protocol mirrors the currently required simulator operations
- the factory can instantiate the current backend path
- no public API changes

## PR 3: C++ Backend Adapter

Priority: P0

Status:

- Complete

Goal:

- hide the current pybind implementation behind the backend protocol

Deliverables:

- `clode/_backends/cpp.py`
- public simulators instantiate the backend via the factory instead of direct pybind construction

Checkpoint:

- the default backend remains the current C++ path

Acceptance criteria:

- `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` still behave the same
- the pinned core suite passes unchanged
- backend-specific implementation details are no longer constructed directly in public modules

## PR 4: RHS Source Object Integration

Priority: P0

Status:

- Complete

Goal:

- make internal backend calls consume source text plus digest rather than only file paths

Deliverables:

- internal `RhsSource` preparation in the public Python layer
- no behavior change for users providing file paths or Python RHS callables

Checkpoint:

- backend seam can pass source text to later PyOpenCL components

Acceptance criteria:

- current C++ backend still works through the adapter
- `get_program_string()` behavior is unchanged on the C++ backend
- converter and XPP tests still pass

Audit result:

- Passed

## PR 5: PyOpenCL Core Models And Errors

Priority: P0

Status:

- Complete

Goal:

- land the core PyOpenCL data models and error types without execution logic

Deliverables:

- `clode/_pyopencl/models.py`
- `clode/_pyopencl/errors.py`
- `clode/_pyopencl/__init__.py`

Checkpoint:

- build-key model is reviewed and accepted before source builder work starts

Acceptance criteria:

- data models cover build key, source bundle, and program bundle needs
- no execution yet
- no public API changes

Audit result:

- Passed

## PR 6: PyOpenCL Registry And Source Builder

Priority: P0

Status:

- Complete

Goal:

- reproduce the current build-key and source assembly logic in Python

Deliverables:

- `clode/_pyopencl/registry.py`
- `clode/_pyopencl/source_builder.py`

Checkpoint:

- transient source assembly works and emits inspectable source plus build options

Acceptance criteria:

- source builder reproduces current transient, trajectory, and feature entrypoint composition
- build options include current compile-time specialization knobs
- RHS validation exists at least at a conservative phase-one level
- tests cover build-key changes and source assembly behavior

Audit result:

- Passed

## PR 7: PyOpenCL Runtime And Program Cache

Priority: P0

Status:

- Complete

Goal:

- compile programs and cache them correctly in Python without executing solver kernels yet

Deliverables:

- `clode/_pyopencl/runtime.py`
- `clode/_pyopencl/program_cache.py`

Checkpoint:

- transient program compiles successfully on a supported device through PyOpenCL

Acceptance criteria:

- runtime selects one device explicitly
- cache is runtime-scoped and keyed by `BuildKey`
- build failures surface source text, options, and build log

Audit result:

- Passed

## PR 8: Buffer Manager For Common State

Priority: P0

Status:

- Complete

Goal:

- centralize common buffer allocation and current flatten and reshape rules

Deliverables:

- `clode/_pyopencl/buffers.py`

Checkpoint:

- common solver buffers exist independently of any specific executor

Acceptance criteria:

- current Fortran-order flattening is centralized
- common buffer allocation supports transient execution requirements
- no public API changes

Audit result:

- Passed

## PR 9: PyOpenCL Transient Backend

Priority: P0

Status:

- Complete

Goal:

- land the first working PyOpenCL execution path

Deliverables:

- `PyOpenCLTransientBackend` in `clode/_pyopencl/executors.py`
- factory support for selecting the PyOpenCL backend internally
- `clode/_pyopencl/structs.py` for device-matched host struct packing used by the transient path
- `.design/pyopencl_struct_audit.md` recording the `ObserverData` sizing risk that must not be copied into later PRs

Checkpoint:

- transient runs locally through PyOpenCL behind an internal backend selector

Acceptance criteria:

- transient parity tests compare C++ backend to PyOpenCL backend
- seeded stochastic repeatability behavior is pinned
- continuation of `x0`, `dt`, and final time matches current behavior

Audit result:

- Passed

## PR 10: Trajectory Buffers And Trajectory Backend

Priority: P1

Status:

- Complete

Goal:

- add trajectory execution without introducing the future monitor abstraction

Deliverables:

- trajectory buffer support in `clode/_pyopencl/buffers.py`
- `PyOpenCLTrajectoryBackend` in `clode/_pyopencl/executors.py`
- `test/test_pyopencl_trajectory_backend.py`

Checkpoint:

- `TrajectoryOutput` remains unchanged

Acceptance criteria:

- trajectory contract tests pass against both backends
- `test/test_ornl_thompson_a1.py` passes on the PyOpenCL path
- `nout`, `max_store`, and `n_stored` semantics match the current backend

Audit result:

- Passed

## PR 11: Feature Buffers And Feature Backend

Priority: P1

Status:

- Complete

Goal:

- add the current observer and feature pipeline to the PyOpenCL backend

Deliverables:

- feature buffer support in `clode/_pyopencl/buffers.py`
- observer metadata and layout sizing in `clode/_pyopencl/observer_metadata.py`
- `PyOpenCLFeatureBackend` in `clode/_pyopencl/executors.py`
- `test/test_pyopencl_feature_backend.py`

Checkpoint:

- one-pass and two-pass observers both execute through PyOpenCL

Acceptance criteria:

- feature contract tests pass against both backends
- `test/test_vdp.py`, `test/test_features.py`, and `test/test_aux_values.py` pass on the PyOpenCL path
- observer initialization and continuation behavior match current behavior

Guardrails:

- do not reuse the legacy C++ `observerDataSize` byte-count formulas for PyOpenCL allocation
- treat `ObserverData` as observer-specific opaque device state whose exact size comes from an explicit layout model
- preserve the current `initializeObserver` plus `features` lifecycle, including the two-pass warmup path used by `nhood2` and `thresh2`

Audit result:

- Passed on the stable NVIDIA runtime; the Intel CPU runtime remains unstable for the `localmax` rebuild path on this workspace

## PR 12: Extended Reference Suite And Rollout Guardrails

Priority: P1

Goal:

- widen validation without changing the default backend yet

Status:

- Complete

Deliverables:

- extended reference suite runs in milestone validation on the stable supported runtime
- rollout docs updated with current backend-selector policy and runtime-selection guardrails
- `.design/pyopencl_rollout_guardrails.md`

Checkpoint:

- `test/test_opencl_builtins.py` and `test/test_runtime.py` are carried into the extended validation bundle

Acceptance criteria:

- the 74-test extended reference bundle is runnable and documented
- backend selector remains internal
- C++ fallback remains intact

Audit result:

- Passed on the stable NVIDIA runtime; the Intel CPU runtime remains outside the current supported rollout set on this workspace

## PR 13: Dependency Surfacing And Runtime Diagnostics

Priority: P1

Status:

- Complete

Goal:

- make the PyOpenCL installation path and runtime-selection story explicit before attempting the default switch

Deliverables:

- package metadata and contributor requirements make the PyOpenCL dependency visible
- install docs describe runtime verification and supported rollout expectations
- runtime diagnostics are clear enough that selected platform and device tuples can be verified without guesswork

Checkpoint:

- maintainers can install the PyOpenCL path intentionally and confirm which runtime was selected

Acceptance criteria:

- the PyOpenCL dependency story is explicit without making the public default backend change yet
- runtime-selection guidance points at `clinfo -l` plus the in-package OpenCL query helpers
- the C++ backend remains the default path during this PR

Audit result:

- Passed on the stable NVIDIA runtime; the 74-test extended reference bundle remained green after the packaging, install-doc, and runtime-diagnostics updates

## PR 14: Python-Owned Public Types And Runtime Facade

Priority: P1

Status:

- Complete

Goal:

- remove the C++ wrapper from the public Python import path so the package can become PyOpenCL-owned rather than only PyOpenCL-capable

Deliverables:

- Python-owned `ProblemInfo`, `SolverParams`, and `ObserverParams`
- Python-owned runtime compatibility models for public runtime operations
- adapter conversions at the C++ backend boundary instead of wrapper-type imports throughout the public package
- `_pyopencl` internals updated to use Python-owned models

Checkpoint:

- `import clode` no longer requires the C++ extension when the PyOpenCL path is installed and selected

Acceptance criteria:

- public modules no longer import `clode.cpp.clode_cpp_wrapper` at module import time
- `_pyopencl` backends no longer depend on C++-owned public structs
- C++ backend still works through explicit adapter conversions

Current slice landed:

- `ProblemInfo`, `SolverParams`, and `ObserverParams` now have Python-owned implementations in `clode/types.py`
- public simulators and `_pyopencl` modules now consume those Python-owned model types
- the C++ backend now converts those models explicitly at the adapter boundary
- `clode/runtime.py` now owns the public runtime enums, info models, lazy logger state, and a Python `OpenCLResource` compatibility facade
- `clode/_backends/factory.py` now lazy-loads the C++ adapter so `import clode` and the PyOpenCL path no longer touch the wrapper at module import time
- `test/test_pyopencl_runtime.py` now contains a subprocess regression that blocks `clode.cpp` imports and verifies the PyOpenCL public path still imports, queries OpenCL, and builds a simulator
- `examples/pyopencl_ornstein_uhlenbeck.py` now demonstrates the current transition-phase PyOpenCL backend selector and explicit runtime pinning on a real simulation
- focused runtime tests, the default-backend runtime/logger smoke, and the broader 76-test PyOpenCL acceptance bundle all passed locally on this workspace

Audit result:

- Passed on the stable NVIDIA runtime; `import clode` no longer requires the C++ extension when the PyOpenCL path is selected, and the C++ backend still works through the explicit adapter boundary

Guardrails:

- do not change kernel numerics or constructor semantics in this PR
- keep the C++ backend working while the public model layer is replaced

## PR 15: Bazel-Free Packaging And Release Transition

Priority: P1

Status:

- completed on this workspace on 2026-05-05

Goal:

- make the default package build and release path pure Python before changing the public backend default

Deliverables:

- default `python -m build` path no longer compiles the C++ extension
- wheel carries Python code plus required OpenCL assets, not the legacy C++ source tree
- stale trees such as `matlab/` and `samples/` are no longer part of package artifacts
- release/version source is unified
- release publishing moves to a tag-driven workflow rather than general push jobs

Checkpoint:

- the default wheel becomes `py3-none-any`

Acceptance criteria:

- `python -m build` succeeds without Bazel on the default package path
- the default wheel is pure Python
- legacy C++ support is retained only through an explicit compatibility path

Validated result:

- `python -m build` now emits `clode-0.9.0-py3-none-any.whl`
- runtime kernel assets are packaged under `clode/kernels/` instead of the legacy `clode/cpp` tree
- the default wheel no longer ships the wrapper extension or the `clode/cpp` source tree
- an isolated wheel install resolves to the PyOpenCL runtime path when the wrapper is absent
- release publishing now lives in one dedicated tag-driven workflow rather than the push CI jobs

Guardrails:

- do not switch the default backend in the same PR
- keep the release workflow understandable enough that published artifacts map cleanly to one semver tag
- use standard `pyproject.toml`-first packaging conventions for the PyOpenCL-only endpoint rather than carrying forward Bazel-era release assumptions

## PR 16: Default Backend Switch

Priority: P1

Status:

- completed on this workspace on 2026-05-06

Goal:

- switch the default backend to PyOpenCL only after the public package and release path are already PyOpenCL-owned

Deliverables:

- backend factory defaults to PyOpenCL
- explicit legacy fallback retained during transition where still supported
- docs updated for runtime and dependency behavior
- source checkouts stop preferring the legacy backend merely because a local wrapper binary is present

Checkpoint:

- maintainers explicitly sign off on parity and packaging readiness

Acceptance criteria:

- the extended reference bundle passes with PyOpenCL as default on supported environments
- default install path no longer depends on the C++ extension build
- no public API changes are introduced in the switch PR

Validated result:

- the backend selector now resolves to PyOpenCL by default even in a source checkout that already contains a locally built wrapper binary
- the legacy C++ backend remains available only through explicit `_CLODE_BACKEND=cpp` selection
- the package version is now `0.10.0`, matching the intended first pure-Python default-backend release line
- public docs now describe the legacy wrapper as a comparison workflow rather than as the default runtime path
- migration-state notes and post-migration cleanup notes remain separated between this file family and `.design/pyopencl_post_migration_plan.md`

Guardrails:

- do not switch the default backend while runtime-selection guidance is still ambiguous on supported machines
- do not switch the default backend until the PyOpenCL installation path is documented strongly enough for new users
- keep an explicit escape hatch to the legacy backend during the transition if it is still shipped

## PR 17 And Later: Post-Parity Simplification

Priority: P2

Goal:

- simplify architecture only after the backend is already shipping

Candidate work items:

1. remove mixed host metadata from `.cl` files
2. unify duplicated kernel prolog and epilog code
3. revisit monitor or sink abstraction
4. separate persistent observer state from optional event-storage capacity
5. add chunked trajectory streaming
6. remove the C++ backend once no longer needed
7. improve runtime-selection diagnostics so `clinfo`, the legacy wrapper, and PyOpenCL are easier to reconcile

Checkpoint:

- every simplification PR must prove reduced complexity, not just moved complexity

Acceptance criteria:

- parity suite remains green
- public contract remains stable unless an explicit follow-on design change is approved

## tmp Doc Status After PR16

Keep active during the remaining migration cleanup:

- `.design/backend_pr_task_breakdown.md`
- `.design/pyopencl_rollout_guardrails.md`
- `.design/pyopencl_packaging_release_audit.md`
- `.design/pyopencl_backend_design.md`

Keep as technical references while the legacy comparison path still exists:

- `.design/cpp_opencl_layer_audit.md`
- `.design/pyopencl_struct_audit.md`
- `.design/python_simulation_flow_reference.md`

Archived on 2026-05-07 during the first wrapper-removal cleanup pass:

- `.design/archived/backend_migration_history_2026_05_07/backend_overhaul_scope_map.md`
- `.design/archived/backend_migration_history_2026_05_07/backend_core_test_suite.md`

Forward-looking notes should continue to live separately in `.design/pyopencl_post_migration_plan.md` and `.design/cpp_wrapper_cleanup_plan.md` rather than being mixed back into the migration-status documents.

## Checkpoint Summary

Use this section for fast progress assessment.

| Checkpoint | Meaning |
| --- | --- |
| C0 | Scope and tests are pinned |
| C1 | Public simulators no longer directly construct pybind backends |
| C2 | Python can compute the build key and assemble source deterministically |
| C3 | Python can compile and cache transient programs through PyOpenCL |
| C4 | Transient parity exists |
| C5 | Trajectory parity exists |
| C6 | Feature parity exists |
| C7 | Extended reference suite passes |
| C8 | Public package no longer requires the C++ wrapper at import time |
| C9 | Default package build is Bazel-free |
| C10 | PyOpenCL is the default backend |

## PR Sizing Rule

If a proposed PR changes more than one of these categories at once, it is probably too large.

- backend seam
- source assembly
- runtime and cache
- common buffers
- transient execution
- trajectory execution
- feature execution
- rollout and default switching

## Maintenance Rule

If the planned order changes, this document must change in the same PR as the first implementation that depends on that new order.
