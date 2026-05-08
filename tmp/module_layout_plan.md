# Module Layout Plan

## Decision

Plan clODE as a PyOpenCL-first package with no planned alternate backends.

That means:

- do not preserve `_backends/` as a long-term architecture seam just to keep hypothetical backend optionality alive
- favor role-based module names in the public package
- keep OpenCL and PyOpenCL implementation details behind underscore-prefixed internal packages
- preserve the current top-level user API during the migration by re-exporting from `clode.__init__`

This is a layout and semantics plan only. It does not require immediate user-facing API changes.

## Design goals

1. Make the public package read like the user mental model, not the migration history.
2. Leave room for the roadmap items that are likely to reshape internals next: continuation state, observer definitions, RHS IR, and output/storage separation.
3. Remove names that imply backend pluralism where the code no longer has it.
4. Keep the migration incremental, with compatibility aliases at the package root.
5. Avoid creating a second round of churn by renaming files without clarifying ownership.

## Constraints from the current roadmap

- Continuation and solver-state work will likely introduce clearer internal state objects.
- Observer cleanup will likely introduce richer observer definitions and cleaner parameter/state separation.
- RHS work will likely grow beyond a single converter file into a more explicit internal representation.
- Runtime selection remains a public concern, but `device_ids` is transitional and real multi-device execution is not planned.
- Kernel assets under `clode/kernels/` should remain the package-data root.

## Current semantic problems in the layout

- The root package is doing too much. `solver.py`, `trajectory.py`, `features.py`, `runtime.py`, `function_converter.py`, `xpp_parser.py`, `types.py`, and `opencl_builtins.py` represent different roles, but they all sit flat at the top level.
- `types.py` is too generic a name for three unrelated configuration and metadata dataclasses.
- `SolverParams` is solver configuration, not generic type infrastructure. It governs how the numerical method advances a problem in time, so it belongs with simulation semantics rather than with problem definition or runtime/device selection.
- `function_converter.py`, `xpp_parser.py`, and `opencl_builtins.py` are really one family: RHS authoring and ingestion.
- `trajectory.py` and `features.py` each mix simulator classes with output container classes.
- `_backends/` still names a choice the package no longer intends to offer.
- `_pyopencl/` names the current binding library rather than the role those modules play in the package.

## Options considered

### Option A: minimal cleanup

Keep the flat public root, keep `_pyopencl/`, and only remove `_backends/`.

Pros:

- smallest migration cost
- least file motion

Cons:

- keeps the flat, mixed-semantics root package
- does not make future state, observer, and RHS work easier to place cleanly
- preserves too much transition history in the layout

Verdict:

- not recommended

### Option B: role-based public subpackages with compatibility re-exports

Introduce public subpackages by role, keep the current top-level imports as compatibility aliases, and collapse `_backends/` into the OpenCL implementation layer.

Pros:

- matches user-facing concepts
- supports future internal growth without more flat-root sprawl
- lets us migrate in phases without breaking the top-level API

Cons:

- more file motion than option A
- requires discipline to avoid exposing too much internal structure too early

Verdict:

- recommended

### Option C: aggressive semantic redesign now

Redesign immediately around explicit problem definitions, simulation state objects, observer definitions, and output-policy objects before moving files.

Pros:

- potentially the cleanest end state

Cons:

- too much simultaneous change while continuation semantics are still unsettled
- likely to entangle layout work with unresolved execution-model work

Verdict:

- good eventual direction for some internals, but too aggressive as the first package-layout move

## Recommended target structure

Keep `clode.__init__` as the stable top-level barrel, but migrate toward this structure underneath it:

```text
clode/
  __init__.py
  problem/
    __init__.py
    definition.py
    python.py
    xpp.py
    builtins.py
    source.py
  runtime/
    __init__.py
    selection.py
    query.py
    logging.py
  simulation/
    __init__.py
    params.py
    base.py
    trajectory.py
    features.py
    results.py
  observers/
    __init__.py
    types.py
    definitions.py
  _opencl/
    __init__.py
    runtime.py
    executors.py
    registry.py
    source_builder.py
    program_cache.py
    buffers.py
    structs.py
    observer_metadata.py
  kernels/
```

## Intended semantics by area

### `problem/`

Public problem-definition and problem-ingestion surface.

Expected contents:

- `definition.py`: `ProblemInfo` now, and a future richer `ProblemDefinition` if that emerges
- `python.py`: Python-authored RHS conversion helpers
- `xpp.py`: XPP parsing and conversion helpers
- `builtins.py`: OpenCL math compatibility names used when authoring Python-side problems
- `source.py`: `RhsSource` and source-loading helpers

Why:

- `ProblemInfo`, source loading, Python-to-OpenCL conversion, and XPP ingestion already meet in `Simulator._prepare_rhs_source`
- this is one user-facing concern: how a model or problem is defined and turned into an OpenCL-ready form
- `rhs/` is an implementer-centric label, while `problem/` matches the public mental model better

Later:

- this is the natural home for any future richer `ProblemDefinition` or stronger RHS/problem IR work

### `runtime/`

Public runtime and device-selection surface.

Expected contents:

- device-selection enums and parsing helpers
- `query_opencl()` and `print_opencl()`
- runtime logging controls

Why:

- users already think of runtime selection and device inspection as one concept
- the current `runtime.py` is large enough that a package is more intuitive than a single file
- `clode.runtime` should be the stable public package, with selection, query, and logging responsibilities split into submodules and re-exported through `runtime/__init__.py`

Notes:

- `OpenCLResource` should stay internal-or-compatibility oriented, not become a promoted public concept

### `simulation/`

Public simulator orchestration surface.

Expected contents:

- `params.py`: `SolverParams`
- `base.py`: `Simulator`, `Stepper`, and future internal state-related helpers
- `trajectory.py`: `TrajectorySimulator`
- `features.py`: `FeatureSimulator`
- `results.py`: `TrajectoryOutput`, `ObserverOutput`

Why:

- these modules all answer the same user question: how to run a model and what kind of result comes back
- `SolverParams` belongs here because it configures the ODE solver itself: time-step bounds, tolerances, storage cadence, and related integration controls are part of how the simulation is carried out
- those controls are distinct from `problem/`, which describes the dynamical system being solved, and from `runtime/`, which describes the OpenCL hardware/runtime used to execute the solver
- `results.py` removes the current split where output-container classes live inside feature- or trajectory-specific files
- `simulation/` should be the canonical implementation home for `Simulator`, `FeatureSimulator`, and `TrajectorySimulator`; the historical root modules `solver.py`, `features.py`, and `trajectory.py` should survive only as compatibility re-exports during the migration
- a small `params.py` module keeps solver configuration adjacent to `Simulator` and `Stepper` without forcing circular imports inside the semantic `simulation/` package

Later:

- this package has a natural home for future explicit solver-state or output-policy abstractions

### `observers/`

Public observer-facing concepts, with room for future richer internal definitions.

Expected contents:

- `types.py`: `Observer`, `ObserverParams`
- `definitions.py`: initially internal-ish support code, later the home for observer-definition cleanup if that becomes public enough to expose

Why:

- observers are a first-class clODE concept, not just a feature-simulator detail
- the roadmap already points toward richer observer definitions and parameter subsets

### `_opencl/`

Internal OpenCL execution layer, implemented through PyOpenCL.

Expected contents:

- runtime ownership, compilation, buffers, struct layout, observer metadata, executors

Why rename `_pyopencl/` to `_opencl`:

- the role is "internal OpenCL execution layer", not "publicly interesting choice of host binding"
- PyOpenCL remains the implementation detail, but no other binding is planned
- the package name should describe its role in clODE, not keep transition history in the path

## What should happen to `_backends/`

Recommended direction:

- remove `_backends/` as a named long-term subsystem
- either inline the tiny remaining seam into `simulation/` plus `_opencl/`, or keep only a very small internal typing helper if it still pulls its weight

Current status:

- `_backends/` no longer acts as a live architectural seam for simulation/runtime code
- canonical execution ownership now lives in `_opencl/`
- `_backends/rhs.py` and the historical factory/protocol paths can remain only as compatibility wrappers until downstream imports and tests stop needing them

Specific consequence:

- `factory.py` should not survive as a major architectural boundary once the layout work begins
- `protocol.py` is only worth keeping if it materially improves typing or tests; otherwise the public simulators can depend on the concrete `_opencl` executor layer
- `rhs.py` belongs with the rest of problem ingestion, not with a backend package

## Naming rules

1. Public modules should be named by role, not by implementation technology.
2. Internal underscore-prefixed packages may mention OpenCL, but should not imply backend plurality unless there is a real supported choice.
3. Avoid generic catch-all names like `types.py` when the contents clearly belong to problem, runtime, simulation, observer, or RHS roles.
4. Keep result container classes out of simulator modules when they are conceptually reusable outputs.
5. Keep the root `clode.__init__` curated and stable even as the real implementation moves underneath it.

## Migration plan

### Phase 1: establish the public structure without breaking imports

- add the new public subpackages and move or re-export symbols into them
- keep the current root-level imports working through compatibility re-exports
- do not change user-facing behavior yet

### Phase 2: move pure-Python public semantics first

- move problem authoring and parsing into `problem/`
- move output container classes into `simulation/results.py`
- move observer-facing enums and params into `observers/`
- move `ProblemInfo` out of `types.py`
- move `SolverParams` out of `types.py` into `simulation/params.py`
- move the actual simulator classes into `simulation/base.py`, `simulation/features.py`, and `simulation/trajectory.py`, leaving the flat root modules as compatibility barrels

Why first:

- these are mostly semantic and layout changes, not execution changes

### Phase 3: split and clean up the runtime surface

- convert the current `runtime.py` into the `runtime/` package
- keep the public API stable through `runtime/__init__.py` and `clode.__init__`
- de-emphasize or retire transitional compatibility concepts such as broad `device_ids`

### Phase 4: collapse transition-era internals

- move `_backends/rhs.py` into `problem/source.py`
- remove or inline `_backends/factory.py`
- move `_pyopencl/` to `_opencl/`

Status after the current migration:

- `problem/`, `simulation/`, `observers/`, and `runtime/` now own the public semantic surface
- `_opencl/` is now the canonical internal execution package
- `_pyopencl/` has been reduced to a compatibility package-level shim

### Phase 5: internal clustering refinements after the first move lands

- only after the first migration settles, decide whether `_opencl/` itself should later be grouped into smaller internal domains such as compile, memory, and observers
- do not over-factor that package in the same PR as the first layout migration

## Non-goals for the first migration

- no immediate public API redesign
- no immediate solver-state redesign
- no immediate observer-definition redesign
- no move of kernel assets out of `clode/kernels/`
- no attempt to make multi-device execution real as part of layout work

## Consensus points to confirm before implementation

1. Use public role-based subpackages rather than keeping a flat root package.
2. Rename `_pyopencl/` to `_opencl/` as the long-term internal execution package.
3. Remove `_backends/` as a named architectural concept during the migration.
4. Keep `clode.__init__` as a stable compatibility barrel during the transition.
5. Treat `problem/`, `simulation/`, `runtime/`, and `observers/` as the main public semantic areas.

## Recommendation

Adopt option B.

It is the cleanest layout that still respects the current package stability constraints. It makes the PyOpenCL-first decision concrete, removes architecture that exists mainly for history, and leaves room for the future roadmap items without forcing those deeper design choices immediately.
