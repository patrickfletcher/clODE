# Backend Overhaul PR Task Breakdown

## Purpose

This document turns the migration design into PR-sized work items with checkpoints and acceptance criteria.

It is intentionally more granular than the phase plan in `tmp/pyopencl_backend_design.md`. The design phases are the roadmap. This document is the execution plan.

## Current Progress

- Completed: PR 0 scope and test lock
- Completed: PR 1 numerical-core and backend-contract test lock
- Completed: PR 2 backend protocol and factory
- Completed: PR 3 C++ backend adapter and public-wrapper rewiring
- Next: PR 4 RHS source object integration

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
| M2 | Land PyOpenCL build primitives | Not started |
| M3 | Land transient parity | Not started |
| M4 | Land trajectory parity | Not started |
| M5 | Land feature parity | Not started |
| M6 | Switch defaults and simplify later | Not started |

## PR 0: Scope And Test Lock

Priority: P0

Goal:

- freeze the authoritative scope and pinned test suite before implementation starts

Deliverables:

- `tmp/backend_overhaul_scope_map.md`
- `tmp/backend_core_test_suite.md`
- `tmp/backend_pr_task_breakdown.md`

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

- new tests are added exactly as scoped in `tmp/backend_core_test_suite.md`
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

- Next

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

## PR 5: PyOpenCL Core Models And Errors

Priority: P0

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

## PR 6: PyOpenCL Registry And Source Builder

Priority: P0

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

## PR 7: PyOpenCL Runtime And Program Cache

Priority: P0

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

## PR 8: Buffer Manager For Common State

Priority: P0

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

## PR 9: PyOpenCL Transient Backend

Priority: P0

Goal:

- land the first working PyOpenCL execution path

Deliverables:

- `PyOpenCLTransientBackend` in `clode/_pyopencl/executors.py`
- factory support for selecting the PyOpenCL backend internally

Checkpoint:

- transient runs locally through PyOpenCL behind an internal backend selector

Acceptance criteria:

- transient parity tests compare C++ backend to PyOpenCL backend
- seeded stochastic repeatability behavior is pinned
- continuation of `x0`, `dt`, and final time matches current behavior

## PR 10: Trajectory Buffers And Trajectory Backend

Priority: P1

Goal:

- add trajectory execution without introducing the future monitor abstraction

Deliverables:

- trajectory buffer support in `clode/_pyopencl/buffers.py`
- `PyOpenCLTrajectoryBackend` in `clode/_pyopencl/executors.py`

Checkpoint:

- `TrajectoryOutput` remains unchanged

Acceptance criteria:

- trajectory contract tests pass against both backends
- `test/test_ornl_thompson_a1.py` passes on the PyOpenCL path
- `nout`, `max_store`, and `n_stored` semantics match the current backend

## PR 11: Feature Buffers And Feature Backend

Priority: P1

Goal:

- add the current observer and feature pipeline to the PyOpenCL backend

Deliverables:

- feature buffer support in `clode/_pyopencl/buffers.py`
- `PyOpenCLFeatureBackend` in `clode/_pyopencl/executors.py`

Checkpoint:

- one-pass and two-pass observers both execute through PyOpenCL

Acceptance criteria:

- feature contract tests pass against both backends
- `test/test_vdp.py`, `test/test_features.py`, and `test/test_aux_values.py` pass on the PyOpenCL path
- observer initialization and continuation behavior match current behavior

## PR 12: Extended Reference Suite And Rollout Guardrails

Priority: P1

Goal:

- widen validation without changing the default backend yet

Deliverables:

- extended reference suite runs in CI or milestone validation
- rollout docs updated with current backend-selector policy

Checkpoint:

- `test/test_opencl_builtins.py` passes through PyOpenCL on supported environments

Acceptance criteria:

- extended suite is runnable and documented
- backend selector remains internal
- C++ fallback remains intact

## PR 13: Default Backend Switch

Priority: P1

Goal:

- switch the default backend to PyOpenCL only after parity is proven

Deliverables:

- backend factory defaults to PyOpenCL
- fallback path to C++ backend retained during transition
- docs updated for runtime and dependency behavior

Checkpoint:

- maintainers explicitly sign off on parity readiness

Acceptance criteria:

- core and extended suites pass with PyOpenCL as default on supported environments
- C++ backend fallback still works
- no public API changes are introduced in the switch PR

## PR 14 And Later: Post-Parity Simplification

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

Checkpoint:

- every simplification PR must prove reduced complexity, not just moved complexity

Acceptance criteria:

- parity suite remains green
- public contract remains stable unless an explicit follow-on design change is approved

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
| C8 | PyOpenCL is the default backend |

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
