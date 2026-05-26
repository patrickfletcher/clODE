# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Solver-owned per-work-item state and failure reporting

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The first-pass Python-side `SolverState`, observer definitions, observer-parameter resolution, stepper definitions, and result-cache invalidation boundaries are live.
- Fixed-step and adaptive live steppers now share one compensated solve-relative elapsed-time story for endpoints and internal stages.
- The wrapper boundary now preserves both stepper status and accepted step width while leaving the next-step proposal in `dt`.
- The runtime and simulator layer now surface solver-owned per-work-item status, accepted-step-count, and last-accepted-step-width arrays, while richer time-base state and any work metrics are still not explicit device-side solver state.
- `features.cl` now feeds observers the accepted step width from the solver boundary rather than recovering it from elapsed-time differences.
- Public observer feature surfaces no longer report step count or `dt` summary diagnostics; remaining observer-private counters only persist where event geometry or internal running means still need them.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.

## Why this should be next

The numerical time-base cleanup has now made the architectural boundary clearer rather than fuzzier: the solver owns solve-relative elapsed time, accepted step width, next-step proposal, and failure policy. The remaining mismatch is that some observer paths and outputs still carry legacy step or time diagnostics that should be solver-owned state instead of observer-owned bookkeeping.

The wrapper-level prep is already in place. The next leverage point is to turn that into an explicit per-work-item solver-state model and a deterministic surfaced failure policy before any deeper observer-memory or continuation redesign. That keeps single source of truth with the solver and lets observers ingest time values rather than invent or report them.

## Scope

- extend the now-landed solver-owned status, accepted-step-count, and last-accepted-step-width path toward any remaining current time-base diagnostics and later work metrics
- make the runtime persist and fetch that solver-owned state separately from observer state
- keep solver-owned diagnostics on solver fetch paths rather than reintroducing them through observer outputs or metadata
- keep only interpolation or event-timestamp sample geometry in observers where it is still semantically observer-owned
- make the failure-to-accept-step policy deterministic and surfaced, even if the public API exposure stays narrow in this PR
- decide whether no-progress from float32 time quantization belongs in the same surfaced failure-status model and, if so, represent it as solver-owned status rather than as a wrapper-side special case
- keep the public API stable unless a small diagnostic accessor is clearly needed and low-risk

## Likely Internal Shape

- choose and document one internal per-work-item solver-state representation for the remaining time-base and work-metric facts, building on the landed status, step-count, and last-accepted-step-width buffers
- thread solver-owned status and stepping outputs through `_opencl/buffers.py`, `_opencl/executors.py`, and the simulation-state layer without widening public config scope
- keep observer public outputs free of solver-owned step or time diagnostics while leaving observer-private counters only where the event logic still needs them
- keep observer buffers focused on event semantics, interpolation geometry, and retained event samples
- defer a bundled public solver-stats object until the remaining fields are stable enough that one object would reduce churn rather than freeze an incomplete story

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, and fixed-stage reconstruction semantics
- do not reintroduce duplicated time or step bookkeeping across solver and observer code paths
- keep any adopted precision-loss or no-progress detection on the same solver-owned status path instead of adding ad hoc wrapper or observer reporting
- keep event count observer-owned and do not reintroduce step-count or `dt` summary outputs on observer public surfaces while the remaining solver-owned diagnostics settle
- do not freeze SciPy-style work metrics or a bundled stats object until their semantics are explicit across fixed, adaptive, stochastic, and any later implicit paths
- keep interpolation buffers and event timestamps in observers where they are still needed for observer semantics
- do not fold trajectory-output policy, observer metadata, and solver-state refactoring into one large rewrite
- do not turn this into a broad performance campaign or a continuation-policy redesign
- keep public docs and design notes explicit about what is solver-owned, what remains observer-owned, and what is still deferred

## Non-goals

- no implicit or IMEX solver work
- no multi-device work
- no broader diverged-time continuation-policy redesign in the same PR
- no broad observer-feature redesign beyond removing solver-owned legacy diagnostics from observers
- no citation metadata, release-tag, or broader repo-surface cleanup in the same PR

## Suggested Implementation Slices

1. Treat the landed status, step-count, and last-accepted-step-width accessors as the stable public slice for this PR.
2. Keep observer public outputs aligned with that boundary while retaining only private counters needed for observer internals.
3. Defer bundled stats and heavier work metrics until the remaining fields have one coherent solver-owned story.

## Code-Facing Checklist

- `clode/simulation/_state.py`, `clode/simulation/base.py`: keep solver-owned execution progress distinct from fetched results and IVP problem data
- `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`: add or formalize per-work-item solver-state buffers and fetch semantics
- `clode/kernels/transient.cl`, `clode/kernels/features.cl`, `clode/kernels/trajectory.cl`, `clode/kernels/initializeObserver.cl`, `clode/kernels/steppers/`: keep solver-owned status, time, and accepted-step-width plumbing consistent
- `clode/kernels/observers/`: remove legacy solver-diagnostic ownership while preserving event and interpolation state
- `test/kernel_components/`, `test/test_simulation_contracts.py`, and the continuation numerics slices: add direct evidence for the new solver-state boundary and failure policy

## Acceptance Criteria

- the runtime has one explicit solver-owned internal home for per-work-item completion or failure status, accepted step counts, and last accepted step width, and any remaining time-base values land on that same path
- observer public outputs no longer own or report solver-state diagnostics except for event count and sample geometry still needed for interpolation or event timestamps
- failure to accept a step has one deterministic policy in the kernels and that status is available to the runtime on fetch
- any adopted no-progress or precision-loss condition uses the same solver-owned status path instead of a parallel wrapper-only signal
- a bundled solver-stats object is intentionally deferred until the remaining fields and work metrics are stable enough to justify freezing one public shape
- the `.design` docs describe the solver/observer boundary in the same way the code now implements it

## Follow-on If This Lands Cleanly

1. numerical validation and evidence bundle refresh, with docs and examples tied to maintained tests
2. observer-state and register-pressure audit once solver-owned diagnostics are no longer mixed into observer state
3. any later per-work-item `t0` or diverged-time continuation follow-through only if that becomes a stronger user-facing need
