# Semantic Layout Audit

Purpose: audit the post-migration package layout for remaining semantic friction, compare against peer solver libraries, and recommend a clearer concept-first direction for future refactors.
Read when: considering module moves, Python-side domain models, contributor-facing layout clarity, or how steppers, solver state, observers, and ODE problems should be represented.
Update when: the recommended layout direction changes, a major semantic model lands, or the external comparison set changes materially.

## Bottom line

- Keep the current top-level semantic packages: `clode.problem`, `clode.observers`, `clode.simulation`, `clode.runtime`, and `clode._opencl`.
- Keep the flat root barrels for now. They remain useful as collaborator signposts while the canonical package layout settles.
- Do not start with a bulk kernel-file move.
- The first-pass `InitialValueProblem`, solver-state, observer-definition, stepper-definition, and output-policy cleanup has landed.
- Keep compile-time build specification separate from runtime state so program-cache keys and rebuild triggers stay explicit rather than leaking through host-side cache invalidation.
- Treat this note as stable layout guidance rather than active sequencing. Use the root `.design` docs for current priorities.
- `_opencl` should stay focused on runtime/build/buffer/dispatch concerns and consume those semantic definitions, rather than continuing to define core concepts through strings, registries, and struct builders.

## Fast path

Stop after this section unless the task needs the full friction inventory or older layout rationale.

- For current package facts, use `.design/package_state.md` first.
- For active delivery scope, use `.design/next_pr.md` and `.design/ideas.md`; this note is not the current task tracker.
- For ordinary package-placement or ownership questions, read `## Bottom line` and `## Recommended direction` and stop there unless a concrete friction point still needs evidence.
- Read `## Current layout audit` only when you need the detailed reasons behind a layout recommendation.
- Historical peer-library comparisons already live in `../archived/reference_cleanup_2026_05_13/semantic_layout_background.md`.

Historical peer-library comparisons and broader layout-option analysis from the earlier, longer version of this note now live in `../archived/reference_cleanup_2026_05_13/semantic_layout_background.md`.

## Current layout audit

### What is already working well

- Public semantic homes are now reasonably clear: `problem` for ingestion/codegen, `observers` for public observer catalog/config, `simulation` for orchestration/results, and `runtime` for user-visible OpenCL selection/query helpers.
- Runtime-specific execution mechanics are concentrated in `_opencl/runtime.py`, `_opencl/source_builder.py`, `_opencl/program_cache.py`, `_opencl/buffers.py`, and `_opencl/executors.py`.
- The backend-migration scaffolding is mostly gone. Public simulators now construct `_opencl` executors directly, and runtime selection is explicitly single-device.
- The kernel assets are packaged cleanly under `clode/kernels/`.

### Semantic friction points still visible in the live code

#### 1. Simulator classes still carry too much semantic weight

- `Simulator` still owns `set_ensemble()`, `set_repeat_ensemble()`, problem-data shaping, cached readbacks, `t_span`, and continuation entry points that forward into `_opencl`.
- `TrajectorySimulator` and `FeatureSimulator` add storage or observer policy, but still inherit the same broad state-owning base object.
- The current problem layer now has an explicit `InitialValueProblem` that owns RHS semantics, defaults, batched inputs, and remembered result shape, but simulator-side compatibility delegates and caches still carry more semantic weight than they should.

What is missing:

- no settled home yet for richer batch generation, broadcasting helpers, and longer-term result-shape policy beyond the current IVP-owned first pass
- the naming distinction is now sharper, but simulator orchestration still carries more cache and continuation-adjacent semantic weight than the long-term split should keep

Consequence:

- the semantic unit of one solved instance is now explicit, but simulator classes still carry too much convenience, cache ownership, and continuation-adjacent state
- later solver-state cleanup still has to peel orchestration concerns away from cached execution state and output state

#### 2. Solver-owned continuation state exists, but mostly as scattered execution state

- Public-facing time-window and cache state still spans `Simulator` fields, executor host arrays, and OpenCL buffers.
- The continuation fixes landed, but the underlying state is still expressed as parallel fields and companion buffers rather than one deliberate state model.

What is missing:

- no fuller first-class continuation-state model beyond the current first-pass `SolverState`
- no explicit per-work-item `t0`
- no clear home for completion or error flags
- no deliberate semantic boundary between solver state and continuation-specific RNG state

Consequence:

- continuation helpers, per-item divergence, batching, and future RNG or implicit-stepper work still have to reason through scattered fields rather than one durable state contract

#### 3. Stepper definition is now clearer, but continuation ergonomics still sit above a low-level state model

- Public stepper names still live in `clode/simulation/base.py` as the `Stepper` enum.
- Built-in stepper traits and build mapping now live in `clode/simulation/_stepper_definitions.py`.
- OpenCL build lookup now routes through that stepper-definition catalog in `clode/_opencl/registry.py` and `clode/_opencl/source_builder.py`.
- Actual algorithm families still live in `clode/kernels/steppers.cl` and `clode/kernels/steppers/*.clh`.

What is still missing:

- no settled continuation-state contract for the case where work-items finish at different attained `tf`
- no stepper-specific public parameter story yet
- no semantic home for future implicit, IMEX, or controller-rich methods beyond the current internal trait model

Consequence:

- the internal source-assembly and validation path is much easier to follow than before
- the next semantic pressure has shifted from stepper mapping to continuation-state semantics and later solver extension work

#### 4. Observer definition is now clearer, but custom/composable observer semantics are still future work

- Public observer names and broad configuration live in `clode/observers/types.py`.
- Public feature naming routes through `clode/observers/metadata.py` on top of the observer-definition catalog.
- Built-in `ObserverDefinition` instances plus `ResolvedObserverSpec` now live in `clode/observers/_definitions.py`.
- Runtime-specific persistent observer-state layout lives in `clode/_opencl/observer_metadata.py`.
- Observer behavior lives in `clode/kernels/observers.cl` and the individual `observer_*.clh` files.

What is still missing:

- no obvious Python-side home for observer functor/event-style logic if the project ever wants more composable or user-definable observers
- no public custom-observer surface or observer-specific public parameter model yet

Consequence:

- the current built-in observer story is much easier to follow than it was before the cleanup
- the remaining semantic pressure has shifted toward stepper definitions and execution-setting resolution rather than observer naming or layout ownership

#### 4a. Output policy should not collapse into the observer metaphor

- It is tempting to describe `trajectory` as a very dense observer that stores one event every `nout` steps, and `transient` as an observer with no retained output.
- That framing is possible in the abstract, but it is semantically awkward for clODE.

What is clearer:

- `transient` is the no-retained-output solve path
- `trajectory` is a retained-sample output policy layered on one solve path
- observers are stateful feature or event-detection contracts that may additionally expose sparse event outputs

Why keep that distinction:

- trajectory storage is defined by output cadence and retained sample layout, not by event semantics
- forcing trajectories into the observer model would blur output-policy concerns with observer-state and feature-schema concerns
- clODE's distinctive observer story is stronger when observers mean online event or feature extraction, not every possible emitted output

#### 5. Core concepts still reach `_opencl` as primitives more often than as richer objects

- observer selection still becomes a resolved name and define pair before `_opencl` consumes it
- compile-time feature choices such as event-storage capacity are still carried as primitive values rather than concept-owned configuration objects

Consequence:

- the runtime layer still does some concept-definition work that would be clearer on the Python semantic side
- this is now more obvious in continuation and output-policy ergonomics than in stepper build mapping

### Things that look transitional but should stay for now

- The flat root barrels. They are still useful as collaborator orientation guides and compatibility shims.
- The simulator classes themselves. The cleanup is about what they compose and own, not about removing them as the public handle.
- The old combined-kernel `odedriver` context. It now lives under `.design/tmp/odedriver.cl` as scratch design material rather than active runtime code.
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

### 6. Keep build specification separate from runtime state

Likely shape:

- an explicit `BuildSpec` or `KernelSpec` carrying only compile-time specialization inputs such as precision, problem shape, stepper family, observer definition, and any remaining compile-time storage-layout choices
- runtime values such as `x0`, parameters, `t_span`, RNG state, and fetched outputs stay out of that build specification
- `_opencl/program_cache.py`, `_opencl/source_builder.py`, and simulator invalidation logic consume that explicit build spec rather than inferring rebuild decisions from mirrored runtime fields

What this would buy:

- clearer program-cache keys and rebuild triggers
- less accidental coupling between runtime invalidation logic and kernel-specialization decisions
- a cleaner separation between semantic state owners and the transfer details that implement them

### 7. Make stepper definitions first-class on the Python side

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

### 8. Build on landed observer definitions rather than reopening observer layout first

Likely shape inside `clode.observers`:

- keep the public enum if useful
- keep `ObserverDefinition` and `ResolvedObserverSpec` as the internal semantic layer for built-in observers
- add a separate authoring story only if custom/composable observers become active work
- keep runtime-specific state layout in `_opencl`, but derive it from the semantic observer definition where possible

What this would buy:

- preserves the landed observer boundary instead of reopening it for another semantics pass
- makes later custom/composable observer work incremental rather than another semantic reset

### 9. Delay kernel relocation until the semantic models exist

This is the key sequencing choice.

The useful first move is not “put the `.clh` files next to some new Python files”. The useful first move is “create Python-side concept definitions that make the kernels easier to understand and assemble”. Once those exist, the project can revisit whether physical co-location of kernel assets would simplify or just reshuffle the include tree.

## Stable implications for future layout work

- Keep the landed solver-state and output-policy boundary as the contract for future layout changes; see `.design/reference/solver_state_implementation_plan.md` for the current boundary record.
- Build on the current execution-setting and observer-definition boundaries rather than reopening the top-level package split.
- Prefer IVP-owned batch helpers, narrower simulator orchestration, and explicit build/runtime separation over inventing new top-level semantic containers prematurely.
- Revisit continuation-state modeling, kernel-math helper sharing, or any future kernel relocation only when a concrete owner-model need justifies it.
- For the current active sequencing and PR target, use `.design/next_pr.md` and `.design/ideas.md`.

## Guidance for future package moves

- Bias new semantic code toward `problem`, `observers`, `simulation`, and a possible `steppers` home.
- Bias runtime/build/dispatch code toward `_opencl`.
- Keep flat barrels as collaborator guides until the semantic layout is more self-evident.
- Avoid adding abstractions that are only justified by hypothetical future backends.
- Prefer moving one concept family at a time over a broad repo reshuffle.
