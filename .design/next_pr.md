# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Selective adoption of evidence-backed time-base and observer helpers

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, observer-parameter, and stepper-definition boundaries are live.
- `Simulator.advance_tspan_to_attained_final_time()` now provides a narrow exact-continuation helper for the shared-final-time case.
- `FeatureSimulator` now resolves legacy `observer_*` compatibility inputs through one canonical `ObserverParams` path and copies caller-provided bundles instead of aliasing them.
- `clode/kernels/clODE_utilities.cl` now carries a shared helper layer with tested prototypes for compensated means, compensated time pairs, fixed-step counter time reconstruction, threshold timestamps, and bounded three-sample extrema.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- `test/kernel_components/` now provides direct OpenCL component coverage for the helper prototypes plus the `basic` and `basicall` observer contracts.
- The current stepper family still derives the next absolute time from the current float32 time plus `dt`, regardless of whether individual kernels use `+=` or an explicit temporary.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The evidence-gathering PR is now strong enough to stop hypothesizing and start adopting the conservative helpers that clearly pay off. The package has direct helper tests, a public accuracy example, and tighter docs; the next step should be a narrow implementation pass, not another broad numerics audit.

This is the right next PR because it converts the clearest wins into live kernel behavior without dragging in the larger adaptive-time or implicit-solver redesigns. It keeps scope narrow, improves real runtime behavior, and leaves the broader dual-realtype adaptive time-base work for a dedicated follow-on.

## Scope

- adopt inverse-linear threshold timestamps in `threshold_2` and any directly analogous observer paths where the crossing is already bracketed by two samples
- adopt the shared bounded three-sample max/min helpers in `local_max` and any directly analogous event bookkeeping that already uses the same three-sample geometry
- implement fixed-step absolute-time reconstruction from `t0 + step * dt` in the fixed-step solver path, using a wider internal step counter where feasible without breaking the current public compatibility surface
- update direct OpenCL component coverage and focused simulator/observer tests for the live helper adoption
- update docs and design notes so they distinguish newly live helper behavior from the still-deferred adaptive-time redesign
- keep the public API stable in this PR
- keep adaptive dual-realtype time-base work, public packaging hygiene, citation metadata, and release-tag cleanup intentionally out of scope

## Likely Internal Shape

- switch the live observer kernels that are ready now from sampled timestamps or sample-pick extrema to the shared conservative helpers
- move the fixed-step solver path off repeated float32 absolute-time addition and onto a counter-reconstructed time path
- update the component tests so the helper prototypes that were added in the evidence PR now also cover the adopted live behavior
- keep the public docs explicit about which strategies are now live and which remain prototype-only or deferred
- leave adaptive dual-realtype time-base implementation, broader per-work-item `t0`, and any helper that would require general nonlinear system solve machinery to later explicit PRs

## Design Constraints

- no device-side per-work-item `t0` redesign in this PR
- no public config redesign in this PR
- preserve the landed solver-state, continuation-helper, output-policy, observer-state, observer-parameter, execution-setting, and stepper-definition boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- do not claim that fixed-step counter reconstruction solves the adaptive-step time-base problem; adaptive steppers still need a richer time representation
- inverse-linear threshold timestamps are the conservative live target; slope-aware Hermite interpolation stays prototype-only until its robustness on ambiguous/noisy crossings is better characterized
- no helper in this PR should depend on a general nonlinear system solver; anything that needs that should be deferred with the later implicit-method work
- public docs should recommend starting autonomous feature windows near `t = 0` when that does not change model semantics
- public docs should distinguish newly live helper adoption from the remaining prototype-only helpers instead of implying broader rollout than the code actually has
- public docs should compare current package features to relevant algorithmic alternatives, not narrate the package's development history
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no broad helper rollout across all remaining observers in the same PR
- no adaptive dual-realtype time-base rollout in the same PR
- no implicit or IMEX solver implementation in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no broader continuation-policy redesign in the same PR
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Adopt inverse-linear threshold timestamps in the live crossing detectors that already have a clean two-sample bracket.
2. Adopt the shared bounded three-sample max/min helpers in `local_max` and the directly analogous max/min bookkeeping paths.
3. Move the fixed-step solver path to counter-reconstructed absolute time and widen the internal step index where feasible.
4. Extend the component and focused simulator tests to pin the adopted live behavior.
5. Update the docs and design notes to reflect the new live helper behavior and the remaining adaptive/implicit deferrals.

## Code-Facing Checklist

- `clode/kernels/clODE_utilities.cl`: keep the shared helper layer as the only implementation home for the adopted interpolation and time-base helpers
- `clode/kernels/observers/observer_threshold_2.clh` and `observer_local_maximum.clh`: adopt only the conservative helpers supported by the current evidence
- `clode/kernels/steppers/` plus the fixed-step driver entrypoints: move fixed-step absolute time off repeated float32 addition
- `test/kernel_components/test_kernel_math.py`: keep direct helper coverage alongside the new live behavior checks
- `docs/numerical_accuracy.md`: explain which helpers are now live and which remain prototype-only
- `docs/examples.md`, `docs/performance_notes.md`, and `docs/continuation.md`: keep the public links and limitation wording aligned
- `.design/reference/single_precision_numerics_note.md` and `.design/reference/continuation_timebase_note.md`: keep the active rationale aligned with the lived implementation

## Acceptance Criteria

- the fixed-step live solver path no longer relies solely on repeated float32 absolute-time addition
- the live threshold observer paths that are changed use the shared inverse-linear timestamp helper rather than sampled endpoint times
- the live local-extremum path that is changed uses the shared bounded three-sample helper rather than duplicating sample-pick logic locally
- component tests and focused simulator tests cover the adopted live behavior in addition to the retained helper prototypes
- the public docs and design notes clearly distinguish the new live helper behavior from the still-deferred adaptive dual-time redesign and prototype-only Hermite work
- the `.design` docs continue to defer broader helper rollout, adaptive dual-realtype time work, broader continuation-policy redesign, implicit-solver/nonlinear-solver work, and packaging cleanup until the numerical/state-model surface settles further

## Follow-on If This Lands Cleanly

If this lands cleanly, the next high-value follow-on should be the adaptive-time companion: a solver-owned dual-realtype compensated time base plus relative elapsed-time bookkeeping for the observers that need it. After that, broader helper rollout can stay selective, and any helper that would need general nonlinear system solve machinery should still wait for the later implicit-method work.
