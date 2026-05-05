# Backend Overhaul Scope Map

## Purpose

This document is the durable low-token lookup map for the PyOpenCL backend overhaul.

It answers three questions quickly:

1. Which files and folders are in scope for the backend overhaul.
2. Which files are authoritative for behavior and architecture decisions.
3. Which files are placeholders, partial implementations, or otherwise non-authoritative and must not drive tests or migration design.

This document should be treated as companion memory for:

- `tmp/cpp_opencl_layer_audit.md`
- `tmp/pyopencl_backend_design.md`
- `tmp/backend_core_test_suite.md`
- `tmp/backend_pr_task_breakdown.md`

## Core Memory Items Required To Land The Overhaul Reliably

These are the minimum facts that must stay stable in project memory throughout the migration.

### Architecture memory

- The C++ layer is the runtime scaffolding to replace, not the numerical core to preserve.
- The OpenCL kernel assets under `clode/cpp` are the behavioral reference for parity.
- Compile-time specialization is required in phase one for precision, stepper, observer, problem dimensions, and stored event count.
- The effective build key includes precision, stepper, observer, `N_VAR`, `N_PAR`, `N_AUX`, `N_WIENER`, `N_STORE_EVENTS`, RHS source digest, and kernel tree digest.
- The current kernel assembly path is entrypoint source concatenation plus OpenCL `#include` resolution from the kernel root.
- The current data layout is Fortran-order flattening on the Python side and variable-major indexing inside kernels.
- The current solver semantics include continuation of `x0`, `dt`, final time, RNG state, and observer state.

### Scope memory

- The public Python API must stay stable through the parity phases.
- The backend overhaul is limited to the runtime and execution layers first.
- The migration must not treat all files in `test/` as authoritative.
- Placeholder and partial files must be explicitly excluded from the pinned migration suite.

### Delivery memory

- The backend seam lands before the PyOpenCL runtime.
- Transient parity lands before trajectory parity.
- Trajectory parity lands before feature parity.
- The broader monitor or sink unification is deferred until after parity.

## Do We Need A Dedicated Scope And Memory Doc

Yes.

The overhaul spans public Python wrappers, C++ glue, OpenCL kernels, build logic, and tests. Without a short scope-memory document, every implementation pass will pay repeated token cost rediscovering:

- where authoritative behavior lives
- which tests are real vs placeholder
- which partial files must be ignored
- which directories are relevant for parity vs later cleanup

This document is that scope-memory layer.

## In-Scope Folders

The following folders are in scope for the backend overhaul.

| Folder | Role in overhaul | Authoritative now | Notes |
| --- | --- | --- | --- |
| `clode/` | Public Python API and future backend seam integration | Yes | First Python integration target |
| `clode/_backends/` | Internal backend seam, factory, and current C++ adapter | Yes | Current transition boundary for backend work |
| `clode/_pyopencl/` | New Python-owned runtime and build implementation | Yes | PR 5 starts the package with model and error foundations |
| `clode/cpp/` | Current backend reference implementation and OpenCL kernel source tree | Yes | Behavioral reference, not long-term runtime target |
| `test/` | Source of the pinned migration suite and excluded-placeholder list | Mixed | Only selected tests are authoritative |
| `tmp/` | Active migration docs, scope docs, test-suite docs, task planning | Yes | Current durable planning memory during the backend prep phase |

## Out-Of-Scope Folders For Phase-One Execution Work

These folders are not in scope for the first phases of the backend overhaul.

| Folder | Reason |
| --- | --- |
| `matlab/` | Separate binding layer, not part of Python backend migration |
| `samples/` | Useful as models and examples, not part of backend runtime implementation |
| `paper/` | Not part of runtime migration |
| `bazel/` | Packaging and build support may matter later, but not for backend seam or parity logic |
| `examples/` | Useful for manual checks, not authoritative for parity |

## In-Scope Files By Responsibility

### Public Python API and orchestration

| File | Responsibility | Why in scope |
| --- | --- | --- |
| `clode/__init__.py` | Public exports | Must remain stable through migration |
| `clode/runtime.py` | Public runtime surface | Will become a compatibility layer over the new runtime |
| `clode/solver.py` | Public simulator orchestration | Now delegates backend construction through the internal factory |
| `clode/trajectory.py` | Public trajectory orchestration and output reshaping | Must preserve existing output semantics through the backend seam |
| `clode/features.py` | Public feature orchestration and observer-facing API | Must preserve existing observer-facing API through the backend seam |
| `clode/function_converter.py` | Python to OpenCL RHS conversion | Must continue to feed the backend source builder |
| `clode/xpp_parser.py` | XPP to OpenCL conversion | Must continue to feed the backend source builder |
| `clode/opencl_builtins.py` | Supported OpenCL builtin wrappers | Important source-generation compatibility surface |

### Internal backend seam

| File | Responsibility | Why in scope |
| --- | --- | --- |
| `clode/_backends/protocol.py` | Internal simulator backend contract | Defines the public-wrapper to backend boundary |
| `clode/_backends/factory.py` | Backend selection and construction | Current transition point for swapping backend implementations |
| `clode/_backends/cpp.py` | C++ adapter behind the backend contract | Current reference implementation path during migration |
| `clode/_backends/rhs.py` | Internal RHS source model and digest helpers | Transition input model for future PyOpenCL source assembly |

### PyOpenCL implementation foundation

| File | Responsibility | Why in scope |
| --- | --- | --- |
| `clode/_pyopencl/models.py` | Immutable build, source, and program models | Foundation for deterministic source builder and runtime cache keys |
| `clode/_pyopencl/errors.py` | Backend-specific validation and build errors | Foundation for parity-grade diagnostics in later PRs |
| `clode/_pyopencl/registry.py` | Static stepper and observer registry plus entrypoint mapping | Phase-one source assembly authority |
| `clode/_pyopencl/source_builder.py` | Deterministic program text and build-option assembly | Phase-one replacement for the current C++ build-input path |
| `clode/_pyopencl/runtime.py` | Explicit PyOpenCL context and queue selection | Phase-one runtime foundation for backend execution |
| `clode/_pyopencl/program_cache.py` | Runtime-scoped compiled program cache | Phase-one compile and cache authority |
| `clode/_pyopencl/buffers.py` | Common-state buffer allocation and layout helpers | Phase-one buffer foundation for transient execution |

### Current C++ backend reference

| File | Responsibility | Why in scope |
| --- | --- | --- |
| `clode/cpp/OpenCLResource.hpp` | Current runtime abstraction | Reference for device/runtime capability behavior |
| `clode/cpp/OpenCLResource.cpp` | Context, queue, program build | Reference for runtime and build behavior |
| `clode/cpp/CLODE.hpp` | Base simulator contract | Reference for base state model |
| `clode/cpp/CLODE.cpp` | Base simulator runtime behavior | Reference for build options, uploads, continuation |
| `clode/cpp/CLODEtrajectory.hpp` | Trajectory specialization contract | Reference for trajectory API and buffer model |
| `clode/cpp/CLODEtrajectory.cpp` | Trajectory specialization behavior | Reference for trajectory allocation and launch behavior |
| `clode/cpp/CLODEfeatures.hpp` | Feature specialization contract | Reference for observer lifecycle and feature API |
| `clode/cpp/CLODEfeatures.cpp` | Feature specialization behavior | Reference for feature allocation, warmup, continuation |
| `clode/cpp/CLODEpython.cpp` | pybind API layer | Reference only until backend seam is in place |
| `clode/cpp/clode_cpp_wrapper.pyi` | Python typing surface for C++ wrapper | Useful for API audit, not the long-term backend shape |

### Current OpenCL kernel source tree

| File or folder | Responsibility | Why in scope |
| --- | --- | --- |
| `clode/cpp/transient.cl` | Base transient kernel entrypoint | Reference behavior for transient parity |
| `clode/cpp/trajectory.cl` | Trajectory kernel entrypoint | Reference behavior for trajectory parity |
| `clode/cpp/initializeObserver.cl` | Observer warmup entrypoint | Reference behavior for two-pass observers |
| `clode/cpp/features.cl` | Feature kernel entrypoint | Reference behavior for feature parity |
| `clode/cpp/realtype.cl` | Precision selection | Required in phase-one build model |
| `clode/cpp/clODE_struct_defs.cl` | Solver parameter struct definitions | Required in phase-one build model |
| `clode/cpp/clODE_utilities.cl` | Shared helper functions | Behavioral reference |
| `clode/cpp/clODE_random.cl` | RNG implementation | Behavioral reference |
| `clode/cpp/steppers.cl` | Stepper include hub and host metadata today | In scope now, host metadata to move later |
| `clode/cpp/steppers/` | Stepper implementations | Must be preserved during parity phases |
| `clode/cpp/observers.cl` | Observer include hub and host metadata today | In scope now, host metadata to move later |
| `clode/cpp/observers/` | Observer implementations | Must be preserved during parity phases |

### Build and packaging files relevant to the migration

| File | Responsibility | Why in scope |
| --- | --- | --- |
| `pyproject.toml` | Python package metadata | Relevant for dependency and editable-install docs |
| `setup.py` | Source build path for current C++ backend | Relevant because editable source installs still use it |
| `pytest.ini` | Test marker definitions | Needed for pinned suite docs |
| `clode/cpp/BUILD` | Current C++ and kernel file grouping | Useful reference, not first-phase implementation target |

## Non-Authoritative Or Excluded Files

These files exist in the repo but must not be treated as authoritative for the backend overhaul.

### Kernel and source placeholders

| File | Status | Why excluded |
| --- | --- | --- |
| `clode/cpp/odedriver.cl` | Partial and currently unused | Explicitly marked unused; not authoritative for current behavior |
| `clode/cpp/observers/observer_trajectory.clh` | Partial concept sketch | Not wired into the active observer registry and not authoritative for current feature behavior |

### Test placeholders, stubs, or non-core tests

| File | Status | Why excluded from the pinned migration suite |
| --- | --- | --- |
| `test/test_trajectory.py` | Empty file | No authoritative behavior |
| `test/test_solver.py` | Placeholder with skipped tests | Explicitly marked not ready |
| `test/test_observers.py` | Debug-only skipped tests | Explicitly marked debug validation only |
| `test/test_vdp_long.py` | Long-running stress test | Useful later, not part of the core gating suite |
| `test/test_clODE_utilities.py` | Notes only, not executable tests | Planning file, not authoritative suite input |
| `test/tests.md` | Test planning notes | Helpful context, not the pinned suite itself |
| `test/test_logger.py` | Logging behavior test | Useful reference, not core backend-overhaul gating |

## Pinned Authoritative Test Files

These files are the current authoritative source pool for backend-phase gating.

### Numerical core gate

- `test/core_numerics/test_transient.py`
- `test/core_numerics/test_trajectory.py`
- `test/core_numerics/test_features_basicall.py`
- `test/core_numerics/test_stochastic.py`

### Backend contract gate

- `test/test_backend_contracts.py`
- `test/test_backend_rhs_source.py`
- `test/test_pyopencl_models.py`
- `test/test_pyopencl_source_builder.py`
- `test/test_pyopencl_runtime.py`
- `test/test_pyopencl_buffers.py`

The exact suite definition is pinned in `tmp/backend_core_test_suite.md`.

## Rapid Lookup Guide

Use this section to answer common implementation questions with minimum token cost.

| Question | Read first |
| --- | --- |
| Where does backend selection happen | `clode/_backends/factory.py` |
| Where is the current backend adapter | `clode/_backends/cpp.py` |
| Where is RHS source preparation modeled | `clode/_backends/rhs.py`, `clode/solver.py` |
| Where are the PyOpenCL build and program models | `clode/_pyopencl/models.py`, `clode/_pyopencl/errors.py` |
| Where are the PyOpenCL registry and source builder | `clode/_pyopencl/registry.py`, `clode/_pyopencl/source_builder.py` |
| Where are the PyOpenCL runtime and compile cache | `clode/_pyopencl/runtime.py`, `clode/_pyopencl/program_cache.py` |
| Where are the PyOpenCL buffer layout and allocation rules | `clode/_pyopencl/buffers.py` |
| What is the current build key | `tmp/cpp_opencl_layer_audit.md`, `clode/cpp/CLODE.cpp` |
| How is the current program assembled | `tmp/cpp_opencl_layer_audit.md`, `clode/cpp/CLODE.cpp`, `clode/cpp/transient.cl`, `clode/cpp/features.cl` |
| Which public APIs must remain stable | `clode/solver.py`, `clode/trajectory.py`, `clode/features.py`, `tmp/pyopencl_backend_design.md` |
| What is the current array layout | `tmp/cpp_opencl_layer_audit.md`, `clode/solver.py`, `clode/trajectory.py`, `clode/features.py` |
| What behavior defines transient parity | `clode/cpp/transient.cl`, `clode/cpp/CLODE.cpp` |
| What behavior defines trajectory parity | `clode/cpp/trajectory.cl`, `clode/cpp/CLODEtrajectory.cpp`, `test/core_numerics/test_trajectory.py` |
| What behavior defines feature parity | `clode/cpp/features.cl`, `clode/cpp/initializeObserver.cl`, `clode/cpp/CLODEfeatures.cpp`, `test/core_numerics/test_features_basicall.py`, `test/test_backend_contracts.py` |
| Which tests must not be used as parity authority | This document and `tmp/backend_core_test_suite.md` |

## Lookup Priority Rules

When context is short, consult files in this order.

1. `tmp/pyopencl_backend_design.md`
2. `tmp/backend_core_test_suite.md`
3. This document
4. `tmp/cpp_opencl_layer_audit.md`
5. The current backend seam under `clode/_backends/`
6. The current reference C++ and OpenCL implementation under `clode/cpp/`

## Maintenance Rule

Whenever a file becomes newly authoritative for the migration, or newly disqualified from the migration suite, this document must be updated in the same PR.
