# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Simulation state and output ownership model

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The first-pass Python-side `SolverState`, observer definitions, observer-parameter resolution, stepper definitions, and result-cache invalidation boundaries are live.
- Fixed-step and adaptive live steppers now share one compensated solve-relative elapsed-time story for endpoints and internal stages.
- The wrapper boundary now preserves both stepper status and accepted step width while leaving the next-step proposal in `dt`.
- The runtime and simulator layer now surface solver-owned per-work-item status, accepted-step-count, and last-accepted-step-width arrays, and now report `NO_PROGRESS` both when a positive requested window collapses to zero in runtime precision and for the maintained in-loop float32 time-stall cases currently covered across transient, trajectory, and feature runs.
- `features.cl` now feeds observers the accepted step width from the solver boundary rather than recovering it from elapsed-time differences.
- Public observer feature surfaces no longer report step count or `dt` summary diagnostics; remaining observer-private counters only persist where event geometry or internal running means still need them.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.
- The wrap-up naming stance is now narrower: `SolverState` means solver-owned live state, `TrajectoryOutput` means fetched retained-sample data, device-side `ObserverState` means persistent observer runtime state, and public `ObserverOutput` means the current feature or event readout emitted by `finalizeFeatures(...)` rather than the persistent observer state object.
- The current test surface now has a maintained public-contract slice for collapsed-window and in-loop `NO_PROGRESS` status behavior plus a first exact stable-linear transient evidence pair: an explicit RK4 global-error convergence slice and an adaptive Dormand-Prince tolerance-refinement slice.
- The release-gating numerics bundle now centers on `test/core_numerics/`, while older top-level workflow and scientific numerics files remain available through a supplemental bundle instead of defining the proof layer.
- `TrajectorySimulator` already keeps trajectory output settings separate at runtime, but `SolverParams` still mixes integration and trajectory-output concerns at the compatibility surface.
- `ObserverParams` already exposes derived runtime and event-output views, but built-in observers, persistent observer state layout, and fetched feature output are still braided together across simulator code and `_opencl` runtime paths.
- `_opencl` already owns buffer allocation, struct packing, and matched observer-state layouts, but it still consumes too many mixed policy bundles and primitive values instead of one explicit Python-owned semantic contract.

## Why this should be next

The hardest remaining blocker is no longer observer footprint. It is that solver state, trajectory output policy, observer runtime settings, event-output capacity, persistent observer state, and fetched outputs are still split across compatibility bundles, simulator subclasses, caches, and `_opencl` transfer or metadata helpers.

That is the right next cut because it clarifies the Python-vs-OpenCL ownership story before more observer work, output ergonomics, or public narrative growth piles on top. It also makes the package's distinctive "observer" workflow easier to explain for a stronger JOSS story and easier to extend with additional built-in or later composable observers.

One important framing decision from the pre-PR audit is now explicit: `transient` and `trajectory` should not be forced under the observer metaphor. Transient is the no-retained-output solve path, trajectory is a retained-sample output policy, and observers remain the stateful event or feature path.

## Scope

- make one explicit Python semantic owner for each of: integration settings, trajectory output policy, observer runtime settings, event-output policy, solver state, persistent observer state, and fetched feature or trajectory outputs
- make `_opencl` consume those models for buffer allocation, struct packing, and invalidation rather than silently define the concepts through mixed bundles or mirrored fields
- keep `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` as orchestration handles while reducing the semantic ownership they still carry today
- preserve the public compatibility inputs for now; this PR should clarify internal ownership first, not ship a broad public API redesign
- leave user-facing custom/composable observers and deeper continuation redesign for follow-on PRs, but shape this slice so those follow-ons are smaller and easier to reason about
- preserve the distinction between retained-output policy and observer semantics rather than trying to model trajectory storage as "just another observer"

## Likely Internal Shape

- split the internal owners under `clode/simulation/params.py`, `clode/observers/types.py`, and `clode/simulation/_state.py` into smaller value-oriented models while keeping `SolverParams` and `ObserverParams` as thin compatibility surfaces
- keep persistent observer state distinct from observer runtime settings and from fetched observer output
- keep trajectory output settings distinct from integration settings and solver state
- keep fetched observer readouts conceptually distinct from persistent `ObserverState`; any direct observer-state fetch API or public readout rename is follow-on work, not part of this wrap-up slice
- align `clode/_opencl/executors.py`, `buffers.py`, `observer_metadata.py`, and `structs.py` to consume those explicit owners instead of mixed compatibility bundles

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, solver-owned status boundary, continuation semantics, and current proof layer
- do not reintroduce solver-owned step, status, or time diagnostics through observer outputs or observer-private state
- keep semantic ownership in Python first and treat `_opencl` as the execution consumer rather than the authoritative definition layer
- no user-facing custom-observer DSL, codegen surface, or broad inheritance redesign in the same PR
- avoid turning this PR into a benchmark or register-pressure campaign; the point is architectural clarity, not immediate throughput claims
- keep public docs and design notes explicit about what is clarified here and what remains deferred

## Non-goals

- no implicit or IMEX solver work
- no multi-device work
- no broader diverged-time continuation-policy redesign or matched device-side current-time model in the same PR
- no public bundled solver-stats object in the same PR
- no broad observer-feature redesign or solver-state refactor beyond the owner split needed for this pass
- no public observer-state fetch API or public rename of `ObserverOutput` in the same PR
- no citation metadata, release-tag, or broader repo-surface cleanup in the same PR

## Suggested Implementation Slices

1. Name and codify the internal owners for settings, solver state, persistent observer state, and fetched outputs.
2. Align feature and trajectory executors plus cache invalidation with those owners.
3. Leave one narrow follow-on target for observer authoring/composition once the owner map is explicit.

## Code-Facing Checklist

- `clode/simulation/params.py`, `clode/observers/types.py`, `clode/simulation/_state.py`: make the internal owner split explicit
- `clode/simulation/features.py`, `clode/simulation/trajectory.py`, `clode/simulation/results.py`: keep simulator orchestration and fetched output handling aligned with those owners
- `clode/_opencl/executors.py`, `clode/_opencl/buffers.py`, `clode/_opencl/observer_metadata.py`, `clode/_opencl/structs.py`: make runtime code consume the semantic models rather than define them first
- `.design/reference/semantic_layout_audit.md`, `.design/reference/solver_state_implementation_plan.md`, and `.design/reference/joss_audit.md`: keep the architecture and publication story aligned with the live owner model

## Acceptance Criteria

- the root `.design` docs and touched code can name one owner for integration settings, trajectory output policy, observer runtime settings, event-output policy, solver state, persistent observer state, and fetched outputs
- the feature and trajectory paths stop relying on mixed compatibility bundles internally when a narrower owner model already exists
- the live docs can state clearly that `ObserverOutput` is currently a readout object distinct from persistent observer state, without committing this PR to a public observer-state fetch surface
- the updated design docs give a coherent answer to what solver state, observer state, and outputs mean on the Python side versus the OpenCL side
- the next follow-on for observer authoring/composition is smaller and more explicit than it is today

## Follow-on If This Lands Cleanly

1. observer authoring and composition follow-through built on the explicit owner split
2. numerical evidence and publication-facing benchmark follow-through for the JOSS story
3. later observer-state/register-pressure work once the semantic model is settled enough that the optimization target is clear
