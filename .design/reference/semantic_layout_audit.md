# Semantic Layout Audit

Purpose: audit the post-migration package layout for remaining semantic friction, compare against peer solver libraries, and recommend a clearer concept-first direction for future refactors.
Read when: considering module moves, Python-side domain models, contributor-facing layout clarity, or how steppers, solver state, observers, and ODE problems should be represented.
Update when: the recommended layout direction changes, a major semantic model lands, or the external comparison set changes materially.

## Bottom line

- Keep the current top-level semantic packages: `clode.problem`, `clode.observers`, `clode.simulation`, `clode.runtime`, and `clode._opencl`.
- Keep the flat root barrels for now. They remain useful as collaborator signposts while the canonical package layout settles.
- Do not start with a bulk kernel-file move.
- The first-pass `InitialValueProblem` cleanup has landed, and the next highest-value semantic cleanup is explicit solver-state ownership so simulators read more clearly as orchestration objects with one state boundary per concern.
- Treat stepper definitions and richer observer definitions as follow-on internal semantic-layer work that should build on that clearer solver-state contract.
- `_opencl` should stay focused on runtime/build/buffer/dispatch concerns and consume those semantic definitions, rather than continuing to define core concepts through strings, registries, and struct builders.

Historical peer-library comparisons and broader layout-option analysis from the earlier, longer version of this note now live in `../archived/reference_cleanup_2026_05_13/semantic_layout_background.md`.

## Current layout audit

### What is already working well

- Public semantic homes are now reasonably clear: `problem` for ingestion/codegen, `observers` for public observer catalog/config, `simulation` for orchestration/results, and `runtime` for user-visible OpenCL selection/query helpers.
- Runtime-specific execution mechanics are concentrated in `_opencl/runtime.py`, `_opencl/source_builder.py`, `_opencl/program_cache.py`, `_opencl/buffers.py`, and `_opencl/executors.py`.
- The backend-migration scaffolding is mostly gone. Public simulators now construct `_opencl` executors directly, and runtime selection is explicitly single-device.
- The kernel assets are packaged cleanly under `clode/kernels/`.

### Semantic friction points still visible in the live code

#### 1. Simulator classes still carry too much semantic weight

- `Simulator` still owns `set_ensemble()`, `set_repeat_ensemble()`, problem-data shaping, cached readbacks, `t_span`, and public continuation helpers that forward into `_opencl`.
- `TrajectorySimulator` and `FeatureSimulator` add storage or observer policy, but still inherit the same broad state-owning base object.
- The current problem layer now has an explicit `InitialValueProblem` that owns RHS semantics, defaults, batched inputs, and remembered result shape, but simulator-side compatibility delegates and caches still carry more semantic weight than they should.

What is missing:

- no settled home yet for richer batch generation, broadcasting helpers, and longer-term result-shape policy beyond the current IVP-owned first pass
- no sharper distinction yet between simulator orchestration, solver execution state, and fetched output state

Consequence:

- the semantic unit of one solved instance is now explicit, but simulator classes still carry too much convenience, cache ownership, and continuation-adjacent state
- later solver-state cleanup still has to peel orchestration concerns away from cached execution state and output state

#### 2. Solver-owned continuation state exists, but mostly as scattered execution state

- Public-facing time-window and cache state still spans `Simulator` fields, executor host arrays, and OpenCL buffers.
- The continuation fixes landed, but the underlying state is still expressed as parallel fields and companion buffers rather than one deliberate state model.

What is missing:

- no first-class Python `SolverState`, `TimeWindow`, or `ContinuationState`
- no explicit per-work-item `t0`
- no clear home for completion or error flags
- no deliberate semantic boundary between solver state and continuation-specific RNG state

Consequence:

- continuation helpers, per-item divergence, batching, and future RNG or implicit-stepper work still have to reason through scattered fields rather than one durable state contract

#### 3. Stepper semantics are split across multiple layers

- Public stepper names live in `clode/simulation/base.py` as the `Stepper` enum.
- OpenCL build mapping lives in `clode/_opencl/registry.py` as string-to-define tables.
- Actual algorithm families live in `clode/kernels/steppers.cl` and `clode/kernels/steppers/*.clh`.

What is missing:

- no Python-side `StepperDefinition` or equivalent concept
- no explicit traits for fixed vs adaptive, explicit vs implicit, deterministic vs stochastic
- no semantic home for future implicit, IMEX, or controller-rich methods

Consequence:

- adding or extending steppers currently means touching a public enum, an internal string registry, and kernel assets without one concept-owning definition object

#### 4. Observer definition and observer state are still spread across public catalog, kernel code, and runtime layout

- Public observer names and broad configuration live in `clode/observers/types.py`.
- Public feature naming lives in `clode/observers/metadata.py`.
- Runtime-specific persistent observer-data layout lives in `clode/_opencl/observer_metadata.py`.
- Observer behavior lives in `clode/kernels/observers.cl` and the individual `observer_*.clh` files.

What is missing:

- no unified `ObserverDefinition` that owns feature schema, warmup requirements, event behavior, and the distinction between persistent state and optional output capacity
- no clear semantic story for `ObserverData` versus a higher-level persistent observer-state concept
- no obvious Python-side home for observer functor/event-style logic if the project ever wants more composable or user-definable observers

Consequence:

- collaborator understanding and extension work still require jumping between public types, metadata helpers, OpenCL dtype builders, and kernel include trees
- the solver now owns time, but observer-state naming and ownership have not yet caught up to that model

#### 5. Core concepts still reach `_opencl` as strings more often than as objects

- stepper and observer selection still become raw strings early and are then interpreted by `_opencl/registry.py` and `SourceBuilder`
- compile-time feature choices such as event-storage capacity are still carried as primitive values rather than concept-owned configuration objects

Consequence:

- the runtime layer is still doing some concept-definition work that would be clearer on the Python semantic side

### Things that look transitional but should stay for now

- The flat root barrels. They are still useful as collaborator orientation guides and compatibility shims.
- The simulator classes themselves. The cleanup is about what they compose and own, not about removing them as the public handle.
- `clode/kernels/odedriver.cl`. It remains reasonable to keep as deferred design context rather than active runtime code.
- The separation between public semantic packages and `_opencl`. That split is still valuable; the issue is not the split itself, but where semantic definitions live.

## Recommended direction

### 1. Keep the current top-level package map

The post-migration top-level split is already close to the right one:

- `problem`: problem definition, ingestion, codegen
- `observers`: observer catalog, definitions, public configuration
- `simulation`: orchestration, execution state, results
- `runtime`: user-visible device query/selection/logging helpers
- `_opencl`: runtime-specific execution, buffers, struct layout, build/cache, dispatch

The main improvement should happen inside those semantic homes, not by inventing a new top-level architecture.

### 2. Keep simulators, but narrow them to orchestration

Likely shape:

- `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` remain the public handles.
- They primarily compose an `InitialValueProblem`, solver state, `_opencl` executor/runtime, and output or observer policy.
- Convenience helpers can remain, but they should delegate to clearer IVP-owned batch helpers and state objects rather than being the semantic home of those concepts.

What this would buy:

- preserves the good public entry point while reducing the amount of semantic weight carried by the simulator classes

### 3. Start with `InitialValueProblem`, not `ProblemDefinition`

Likely shape:

- one IVP object owns the RHS semantics plus one mapping of default parameter values and one mapping of default initial-state values
- `ProblemInfo`, `ProblemShape`, and `RhsSource` remain useful as derived static metadata and source/build helpers
- if a lower-level shared definition object later becomes useful, it should probably stay derived or internal rather than being the first new public concept

What this would buy:

- matches the current API, which already requires `variables: Dict[str, float]` and `parameters: Dict[str, float]`
- makes the semantic unit of one solved work item explicit
- gives batch helpers and remembered result shape a clearer home without forcing a second public container type immediately

### 4. Start with IVP-owned batch semantics before splitting out a separate ensemble type

Likely shape:

- let `InitialValueProblem` cover the ordinary size-`(1,)` case plus basic batched parameter or initial-state inputs
- keep shape metadata with that IVP or a very small adjacent value layer so result reshaping remains easy for workflows like parameter grids
- add helper functions for IVP batching, shaping, broadcasting, and generation before deciding whether a dedicated `IVPEnsemble` type adds enough value

What this would buy:

- better semantic clarity around one of clODE's core use cases without committing too early to an extra public container type
- stays aligned with OpenCL build specialization by shared RHS/schema and problem shape
- a cleaner path for grids, random sampling, quasi-random generation, and later batching

When to split later:

- if batch-specific behavior grows beyond normalized arrays plus shape metadata
- if slicing, named batch axes, metadata-rich grid generation, or batching policies become substantial enough that they would make the IVP concept harder to scan

### 5. Model per-work-item solver state separately from solver configuration and results

Likely shape inside `clode.simulation`:

- a small state model around requested window, attained state, per-work-item `t0`, `tf`, `dt`, status flags, cached outputs, and maybe execution status
- a deliberate companion for continuation-specific RNG state if that proves clearer than stuffing everything into one object
- continue using `_opencl/executors.py` and `_opencl/buffers.py` for runtime-specific host/device bookkeeping, but make them implement a clearer semantic state contract

What this would buy:

- solver-owned time base becomes explicit
- continuation helpers and per-item divergence become easier to reason about
- kernel signatures can move toward clearer state-oriented inputs rather than a growing list of parallel buffers

### 6. Make stepper definitions first-class on the Python side

Likely shape:

- `Stepper` enum stays or moves to a stepper-focused module
- add `StepperDefinition` or equivalent with explicit traits such as:
  - fixed vs adaptive
  - explicit vs implicit
  - deterministic vs stochastic
  - kernel define name / kernel family / maybe controller support

What this would buy:

- cleaner implicit-stepper groundwork
- better source assembly than raw string registries
- clearer documentation and contributor mental model

### 7. Promote observers from “enum plus params plus kernel special cases” to definitions

Likely shape inside `clode.observers`:

- keep the public enum if useful
- add a definition object that owns:
  - feature schema
  - warmup/two-pass requirement
  - event semantics
  - public configuration schema
- decide whether a higher-level `ObserverState` concept should sit above the current kernel-side `ObserverData` naming
- keep runtime-specific state layout in `_opencl`, but derive it from the semantic observer definition where possible

What this would buy:

- easier custom/composable observer work later
- better separation between persistent observer state and optional event-output capacity

### 8. Delay kernel relocation until the semantic models exist

This is the key sequencing choice.

The useful first move is not “put the `.clh` files next to some new Python files”. The useful first move is “create Python-side concept definitions that make the kernels easier to understand and assemble”. Once those exist, the project can revisit whether physical co-location of kernel assets would simplify or just reshuffle the include tree.

## Productive near-term follow-on order

1. Make lower-level solver state explicit so continuation helpers and diverged-work-item behavior have a clearer contract.
2. Separate integration state from output and storage policy once the state model has a clearer home.
3. Introduce a Python-owned stepper-definition model and refactor observer definitions around clearer semantic objects.
4. Revisit whether IVP-owned batch helpers are enough or whether a dedicated ensemble type adds real value.
5. Revisit kernel relocation only after those semantic models exist.

## Guidance for future package moves

- Bias new semantic code toward `problem`, `observers`, `simulation`, and a possible `steppers` home.
- Bias runtime/build/dispatch code toward `_opencl`.
- Keep flat barrels as collaborator guides until the semantic layout is more self-evident.
- Avoid adding abstractions that are only justified by hypothetical future backends.
- Prefer moving one concept family at a time over a broad repo reshuffle.
