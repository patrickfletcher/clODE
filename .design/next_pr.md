# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Empirical single-precision numerics demonstrations and guidance

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, observer-parameter, and stepper-definition boundaries are live.
- `Simulator.advance_tspan_to_attained_final_time()` now provides a narrow exact-continuation helper for the shared-final-time case.
- `FeatureSimulator` now resolves legacy `observer_*` compatibility inputs through one canonical `ObserverParams` path and copies caller-provided bundles instead of aliasing them.
- `clode/kernels/clODE_utilities.cl` now carries a small shared compensated-accumulation helper layer.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- `test/kernel_components/` now provides direct OpenCL component coverage for helper math plus `basic` and `basicall` observer contracts.
- The current stepper family still derives the next absolute time from the current float32 time plus `dt`, regardless of whether individual kernels use `+=` or an explicit temporary.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The shared helper layer is now real rather than aspirational, but broader rollout should be evidence-driven. The package still lacks public, reproducible demonstrations that show when compensated mean accumulation materially beats the old `runningMeanTime(...)` path and what candidate time-update strategies do or do not fix.

This is the right next PR because it turns those mitigation ideas into something scientifically inspectable before another broad kernel pass. It improves the package's explanatory value, sharpens later observer and stepper choices, and reduces the risk of overclaiming what the current time-base work actually solves.

## Scope

- add reproducible examples that compare `runningMeanTime(...)` against compensated integral accumulation on float32 workloads that expose late-contribution loss
- add reproducible examples that compare direct time accumulation against compensated or structured time updates and show the limits imposed by float32 spacing, including the limits of reconstructing time from a step counter and the separate large-origin elapsed-time issue for feature windows
- add reproducible observer-safeguard examples for oscillation-oriented workflows, especially amplitude floors, Schmitt-trigger-style threshold hysteresis, and derivative thresholds where they are demonstrably useful
- add reproducible threshold-crossing timestamp comparisons for sampled, inverse-linear, and slope-aware interpolation options
- add a small local-extremum prototype that tests whether the current three-sample buffer geometry supports a meaningfully better interpolation helper
- add narrow kernel-side helper prototypes and component tests for fixed-step counter time, compensated time pairs, threshold interpolation, and bounded three-sample extrema without broad live-kernel rollout in this PR
- add or extend public docs so the relevant float32 theory, current mitigations, and remaining limitations are easy to inspect
- keep the public API stable in this PR
- keep public packaging hygiene, citation metadata, and release-tag cleanup intentionally out of scope

## Likely Internal Shape

- add a focused example script under `examples/` that reproduces the mean and time accumulation issues with float32 arithmetic using the same formulas clODE uses in-kernel
- add a public docs page that explains `ulp(t)`, late-update loss in running means, feature-window origin pitfalls, threshold/local-extremum timestamp tradeoffs, and the current mitigation tradeoffs without overstating them
- add direct OpenCL component tests that pin the new helper prototypes before any later observer or stepper adoption pass
- update the design notes so the later helper-adoption pass is explicitly downstream of this evidence-gathering work
- leave broader per-work-item `t0`, dual-realtype time-base implementation, full structured-time redesign, and any broad kernel-helper rollout to later explicit PRs

## Design Constraints

- no device-side per-work-item `t0` redesign in this PR
- no public config redesign in this PR
- preserve the landed solver-state, continuation-helper, output-policy, observer-state, observer-parameter, execution-setting, and stepper-definition boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- do not claim that `ti = ti + dt` is a mitigation distinct from `ti += dt`; both are the same float32 addition issue
- do not claim that compensated or structured time updates fully solve sub-`ulp(t)` absolute-time resolution in float32
- public docs should recommend starting autonomous feature windows near `t = 0` when that does not change model semantics
- public docs should distinguish current sampled observer timestamps from candidate interpolation helpers instead of implying those helpers are already live everywhere
- public docs should compare current package features to relevant algorithmic alternatives, not narrate the package's development history
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no broad helper rollout across all remaining observers in the same PR
- no implicit or IMEX solver implementation in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no broader continuation-policy redesign in the same PR
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Build a reproducible float32 demonstration for time-weighted means that compares the old recurrence against compensated integral accumulation.
2. Build a reproducible float32 demonstration for time accumulation that compares direct addition against compensated or structured updates, including `step * dt` / `t0 + step * dt` variants.
3. Build a reproducible threshold/local-extremum interpolation comparison and mirror the same ideas in direct OpenCL helper tests.
4. Explain the relevant float32 theory and the observed tradeoffs in public docs.
5. Record the deferred follow-ons explicitly: broader helper adoption, deeper structured-time work, step-counter width and overflow budgets, threshold/local-extremum helper adoption, rounding-mode questions, per-work-item `t0`, adaptive-controller work, and later repo-surface packaging hygiene.

## Code-Facing Checklist

- `examples/single_precision_accuracy.py`: keep the demonstrations small, reproducible, and tied to the exact float32 formulas under discussion
- `docs/numerical_accuracy.md`: explain where the current mitigations help and where float32 representability still wins
- `test/kernel_components/test_kernel_math.py`: pin helper-layer prototypes directly before broad observer or stepper adoption
- `docs/examples.md`, `docs/performance_notes.md`, and `docs/continuation.md`: keep the public links and limitation wording aligned
- `.design/reference/single_precision_numerics_note.md` and `.design/reference/continuation_timebase_note.md`: keep the active rationale aligned with the lived implementation

## Acceptance Criteria

- the repository includes at least one reproducible mean-accumulation demonstration, one reproducible time-accumulation demonstration, one reproducible large-origin elapsed-time demonstration for feature windows, one reproducible threshold-timestamp interpolation demonstration, and one reproducible observer-safeguard demonstration that exercise the current float32 formulas or observer logic directly
- the public docs explain why compensated integral accumulation helps where a naive incremental mean can fail, and they explain why the stepper time-base issue is about float32 representability rather than notation and why large absolute feature-window origins are risky in float32
- the helper layer includes direct component-tested prototypes for fixed-step counter time, compensated time pairs, threshold interpolation, and three-sample extremum recovery, but the PR still avoids claiming broad live-kernel adoption where that work has not happened yet
- the public docs and design notes clearly state what compensated or structured time updates do and do not fix
- the public docs compare current package behavior to relevant algorithmic alternatives without slipping into maintainer-history framing
- the `.design` docs continue to defer broader helper rollout, deeper structured-time redesign, broader continuation-policy work, and packaging cleanup until the empirical guidance is in place and the numerical/state-model surface settles further

## Follow-on If This Lands Cleanly

If this lands cleanly, the next high-value follow-on should be a narrower helper-adoption PR that uses these demonstrations and helper tests to justify where compensated accumulation or time-base changes belong. The preferred single-precision direction after that is fixed-step counter reconstruction where it applies, a solver-owned dual-realtype compensated time base plus relative elapsed-time bookkeeping where it does not, inverse-linear threshold timestamps as the conservative first observer upgrade, and bounded three-sample extremum helpers where the observer contracts benefit clearly.
