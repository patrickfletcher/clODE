# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Numerical validation and evidence bundle

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The first-pass Python-side `SolverState`, observer definitions, observer-parameter resolution, stepper definitions, and result-cache invalidation boundaries are live.
- Fixed-step and adaptive live steppers now share one compensated solve-relative elapsed-time story for endpoints and internal stages.
- The wrapper boundary now preserves both stepper status and accepted step width while leaving the next-step proposal in `dt`.
- The runtime and simulator layer now surface solver-owned per-work-item status, accepted-step-count, and last-accepted-step-width arrays, and now report `NO_PROGRESS` both when a positive requested window collapses to zero in runtime precision and for the maintained in-loop float32 time-stall cases currently covered across transient, trajectory, and feature runs.
- `features.cl` now feeds observers the accepted step width from the solver boundary rather than recovering it from elapsed-time differences.
- Public observer feature surfaces no longer report step count or `dt` summary diagnostics; remaining observer-private counters only persist where event geometry or internal running means still need them.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.
- The current test surface now has a maintained public-contract slice for collapsed-window and in-loop `NO_PROGRESS` status behavior, but the broader numerical evidence layer is still thinner than the current public claims.

## Why this should be next

The next highest-value gap is now proof rather than plumbing. The solver-owned status boundary and the currently adopted `NO_PROGRESS` policy are maintained, but the repo still lacks a slim, explicit solver-validation bundle with clear global-error or convergence expectations that matches the current numerical claims.

That is the right next cut because it raises confidence in the current package story without reopening settled ownership boundaries or mixing in another runtime refactor. It also gives later docs, examples, and performance notes maintained evidence to point at instead of ad hoc demonstrations.

## Scope

- add a slim exact-solution solver-validation slice in `test/core_numerics/` with explicit global-error or convergence expectations for a small representative problem set
- keep release-gating evidence separate from non-gating work-precision or exploratory numerical demos
- align any touched docs or examples with the maintained evidence instead of broadening public numerical claims beyond what the tests prove
- preserve the current solver-owned diagnostics and continuation semantics while validating them through narrower numerical contracts rather than new public API
- keep the public API stable and avoid turning this PR into a broad benchmark or performance campaign

## Likely Internal Shape

- choose one or two realistic ODE problems with exact or high-confidence references and make their acceptance criteria explicit
- keep the authoritative correctness layer in `test/core_numerics/` and use docs or examples only as supplementary demonstrations
- add the smallest helper or reference scaffolding needed for maintained global-error and convergence checks
- leave broader performance plots, benchmark trees, and heavier work-precision experiments outside the release gate

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, solver-owned status boundary, and continuation semantics
- keep correctness evidence separate from performance exploration and avoid turning the test suite into a benchmark harness
- prefer maintained exact-solution or explicit-reference checks over qualitative demo-only evidence
- keep docs and examples explicit about what the maintained tests prove and what remains illustrative only
- do not use this PR to reopen observer-state ownership, current-time modeling, or public solver-stats design

## Non-goals

- no implicit or IMEX solver work
- no multi-device work
- no broader diverged-time continuation-policy redesign or matched device-side current-time model in the same PR
- no broad observer-feature redesign or solver-state refactor in the same PR
- no citation metadata, release-tag, or broader repo-surface cleanup in the same PR

## Suggested Implementation Slices

1. Pick one representative exact-solution or high-confidence reference problem already close to the live workflows.
2. Add explicit global-error or convergence assertions in `test/core_numerics/`.
3. Sync any touched numerical docs or examples so they point to maintained evidence rather than broader claims.

## Code-Facing Checklist

- `test/core_numerics/`: add one narrow solver-validation slice with explicit expectations
- `test/test_ornl_thompson_a1.py`, `examples/`, and `docs/numerical_accuracy.md`: reuse or realign existing evidence where it reduces duplication
- `.design/reference/testing_audit.md` and the root `.design` docs: keep the evidence story aligned with what the maintained tests now prove

## Acceptance Criteria

- at least one representative solver-validation problem has explicit maintained error expectations in `test/core_numerics/`
- the maintained numerical evidence is clearly separated from illustrative examples or work-precision experiments
- the `.design` docs and any touched public numerical docs describe the same evidence boundary the code now enforces

## Follow-on If This Lands Cleanly

1. observer-state and register-pressure audit once the proof layer is tighter
2. examples or docs that show where compensated safeguards help in practice without widening the release gate
3. any later per-work-item current-time or `t0` follow-through only if exact diverged-time continuation becomes a stronger user-facing need
