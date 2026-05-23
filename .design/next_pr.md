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
- `features.cl` now feeds observers the accepted step width from the solver boundary rather than recovering it from elapsed-time differences.
- Observer kernels still retain some legacy step or time diagnostics even though those values semantically belong to solver state rather than to observers.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.

## Why this should be next

The numerical time-base cleanup has now made the architectural boundary clearer rather than fuzzier: the solver owns solve-relative elapsed time, accepted step width, next-step proposal, and failure policy. The remaining mismatch is that some observer paths and outputs still carry legacy step or time diagnostics that should be solver-owned state instead of observer-owned bookkeeping.

The wrapper-level prep is already in place. The next leverage point is to turn that into an explicit per-work-item solver-state model and a deterministic surfaced failure policy before any deeper observer-memory or continuation redesign. That keeps single source of truth with the solver and lets observers ingest time values rather than invent or report them.

## Scope

- define one explicit internal per-work-item solver-state layout for step counts, current or accepted `dt`, current time-base values, and completion or failure status
- make the runtime persist and fetch that solver-owned state separately from observer state
- move legacy step or time diagnostic ownership out of observer structs and outputs where the value is really solver state
- keep only interpolation or event-timestamp sample geometry in observers where it is still semantically observer-owned
- make the failure-to-accept-step policy deterministic and surfaced, even if the public API exposure stays narrow in this PR
- keep the public API stable unless a small diagnostic accessor is clearly needed and low-risk

## Likely Internal Shape

- choose and document one internal per-work-item solver-state representation, likely struct-of-arrays or an equivalent explicit buffer bundle, instead of scattering step or time facts across observers and executors
- thread solver-owned status and stepping outputs through `_opencl/buffers.py`, `_opencl/executors.py`, and the simulation-state layer without widening public config scope
- remove or retire observer-owned step or time diagnostics that are no longer semantically observer state
- keep observer buffers focused on event semantics, interpolation geometry, and retained event samples

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, and fixed-stage reconstruction semantics
- do not reintroduce duplicated time or step bookkeeping across solver and observer code paths
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

1. Define the internal solver-state fields and runtime ownership boundary.
2. Thread solver-owned step, time, and status data through executors and simulation state.
3. Remove or deprecate the matching legacy observer-owned diagnostics.
4. Add direct tests for failure status, step counters, accepted-step-width handoff, and fetched solver-state semantics.

## Code-Facing Checklist

- `clode/simulation/_state.py`, `clode/simulation/base.py`: keep solver-owned execution progress distinct from fetched results and IVP problem data
- `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`: add or formalize per-work-item solver-state buffers and fetch semantics
- `clode/kernels/transient.cl`, `clode/kernels/features.cl`, `clode/kernels/trajectory.cl`, `clode/kernels/initializeObserver.cl`, `clode/kernels/steppers/`: keep solver-owned status, time, and accepted-step-width plumbing consistent
- `clode/kernels/observers/`: remove legacy solver-diagnostic ownership while preserving event and interpolation state
- `test/kernel_components/`, `test/test_simulation_contracts.py`, and the continuation numerics slices: add direct evidence for the new solver-state boundary and failure policy

## Acceptance Criteria

- the runtime has one explicit solver-owned internal home for per-work-item step counts, accepted or current `dt`, time-base values, and completion or failure status
- observers no longer own or report solver-state diagnostics except for sample geometry still needed for interpolation or event timestamps
- failure to accept a step has one deterministic policy in the kernels and that status is available to the runtime on fetch
- the `.design` docs describe the solver/observer boundary in the same way the code now implements it

## Follow-on If This Lands Cleanly

1. numerical validation and evidence bundle refresh, with docs and examples tied to maintained tests
2. observer-state and register-pressure audit once solver-owned diagnostics are no longer mixed into observer state
3. any later per-work-item `t0` or diverged-time continuation follow-through only if that becomes a stronger user-facing need
