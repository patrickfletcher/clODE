# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Numerical validation and evidence bundle

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, observer-parameter, stepper-definition, and result-cache invalidation boundaries are live.
- `Simulator.advance_tspan_to_attained_final_time()` now provides a narrow exact-continuation helper for the shared-final-time case.
- `FeatureSimulator` now resolves legacy `observer_*` compatibility inputs through one canonical `ObserverParams` path and copies caller-provided bundles instead of aliasing them.
- `clode/kernels/clODE_utilities.cl` now carries the shared helper layer for compensated means, compensated time pairs, fixed-step counter time reconstruction, threshold timestamps, and bounded three-sample extrema.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- Fixed-step live steppers now reconstruct absolute time from `t0 + step * dt` with a 64-bit step counter, and adaptive steppers now keep solve-relative compensated elapsed time.
- `threshold_2` now uses inverse-linear timestamps for stored up/down threshold transitions, `local_max` now uses the shared bounded three-sample max/min helpers, and the observer metadata layer now matches the live kernel structs.
- Direct tests now cover program-cache build keys, `t_span` invalidation, observer-layout rebuilds, event-storage changes, and the main helper contracts already adopted in the kernels.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The code is now ahead of the planning docs: the recent execution-state and invalidation cleanup is already landed and covered by tests. The next leverage point is no longer another cache-ownership refactor. It is tightening the proof layer so the numerics, continuation behavior, and public docs say only what the repository can currently demonstrate.

The recent time-base and observer work improved correctness, but it also widened the set of claims around single-precision behavior, split-window continuation, and large-ensemble feature extraction. A narrow validation-and-evidence PR gives those claims a stable backbone before deeper observer-memory or performance architecture work.

## Scope

- audit the current numerical and continuation claims against the live tests, examples, and docs
- add the smallest missing exact-solution, convergence, or component regressions needed to support the landed helper and observer work
- update public docs and design notes so each important claim maps to runnable evidence and the remaining limitations stay explicit
- allow only small, evidence-backed hot-path tidy-ups while keeping the public API stable
- keep the public API stable in this PR
- keep broader observer-time redesign, large-ensemble batching, public packaging hygiene, citation metadata, and release-tag cleanup intentionally out of scope

## Likely Internal Shape

- extend `test/core_numerics/` and `test/kernel_components/` only where the current public or design wording still depends on thin evidence
- keep the docs tied to runnable scripts and direct regressions rather than anecdotal claims or internal narratives
- only touch observer or kernel code when it removes obvious ambiguity or avoidable overhead without expanding the architectural scope
- keep the public docs explicit about what continuation, timestamps, and float32 safeguards do today without widening the public configuration surface in the same PR

## Design Constraints

- no public config redesign in this PR
- preserve the landed solver-state, continuation-helper, output-policy, observer-state, observer-parameter, execution-setting, stepper-definition, and adaptive-time boundaries
- do not turn this into a broad performance campaign or a benchmark-paper pass
- do not reopen execution-state ownership unless a validation gap exposes a real bug
- keep deeper observer-memory or register-pressure redesign on the ideas board unless a small change has direct evidence and low risk
- no helper in this PR should depend on a general nonlinear system solver; anything that needs that should still be deferred with the later implicit-method work
- public docs should stay clear about what continues automatically, what time-base safeguards are live, and where caller-managed policy still exists
- public docs should compare current package features to relevant algorithmic alternatives, not narrate the package's development history
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no public solver-state object redesign in the same PR
- no implicit or IMEX solver implementation in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no broader diverged-time continuation-policy redesign in the same PR
- no deeper observer time-window architecture redesign in the same PR
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Audit the current numerical and continuation claims against the live tests, examples, and docs.
2. Add or tighten the smallest missing exact-solution, convergence, or component regressions.
3. Refresh the public docs and examples so each important claim points to live evidence.
4. Record deeper observer-memory or architecture ideas in `.design/ideas.md` instead of widening the PR.

## Code-Facing Checklist

- `test/core_numerics/` and `test/kernel_components/`: add or tighten the smallest missing direct evidence around the live helper and observer contracts
- `docs/numerical_accuracy.md`, `docs/continuation.md`, `docs/feature_extraction.md`, and `docs/performance_notes.md`: keep user-facing claims aligned with live behavior and reproducible scripts
- `.design/package_state.md`, `.design/ideas.md`, `.design/development_roadmap.md`, and `.design/reference/single_precision_numerics_note.md`: keep the maintainer rationale aligned with the lived implementation and the deferred observer-footprint ideas
- `clode/kernels/observers/`: only take low-risk tidy-ups that clearly reduce avoidable scratch or ambiguity

## Acceptance Criteria

- the main numerical and continuation claims in the docs map to runnable examples or direct tests
- the current helper and observer mitigations have a slim maintained evidence layer rather than only ad hoc regressions
- the public docs and design notes clearly describe what is robust today and what is still limited by float32 spacing or deferred continuation work
- the `.design` docs no longer describe already-landed execution-state or invalidation cleanup as the active next PR
- deeper observer-footprint redesign remains documented as future work rather than being silently expanded into this PR

## Follow-on If This Lands Cleanly

If this lands cleanly, the next grouped follow-ons should stay selective and ordered:

1. observer-state and register-pressure audit, including whether the heavier event detectors can carry a leaner time-window representation without losing the large-origin safeguards
2. large-ensemble ergonomics: IVP-side batch-generation helpers and device-capacity ensemble batching before broader trajectory-output expansion
3. any later per-work-item current-time or `t0` follow-through only if diverged-time continuation becomes a real user-facing priority

Keep broader observer-helper rollout, source-assembly reshaping, broader PyOpenCL helper leverage, and any helper that would need general nonlinear system solve machinery deferred until those grouped follow-ons settle.
