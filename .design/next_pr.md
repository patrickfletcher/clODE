# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Continuation-state semantics and divergence guardrails

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, observer-parameter, and stepper-definition boundaries are live.
- `FeatureSimulator` now resolves legacy `observer_*` compatibility inputs through one canonical `ObserverParams` path and copies caller-provided bundles instead of aliasing them.
- `clode/kernels/clODE_utilities.cl` now carries a small shared compensated-accumulation helper layer.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- `test/kernel_components/` now provides direct OpenCL component coverage for helper math plus `basic` and `basicall` observer contracts.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The main remaining semantic mismatch is now the continuation time model. Each work item owns an attained `tf`, but the current runtime still exposes one shared requested `t_span`. That means exact continuation after diverged per-work-item final times is not representable as one next shared window.

This is the right next PR because it closes a real correctness and UX gap without forcing a premature public continuation API. It should make the internal contract explicit before any richer continuation helpers, batching policy, or trajectory-output redesign build on top of it.

## Scope

- codify that exact shared-window continuation is only valid when the ensemble agrees on one attained `tf`
- add small internal guardrails, queries, or contract tests around that representability limit instead of inventing a broad new public continuation surface
- keep the current low-level `set_tspan(...)`, `shift_tspan()`, and `get_final_time()` API stable while making the unsupported divergent-time case explicit in tests and docs
- keep public packaging hygiene, citation metadata, and release-tag cleanup intentionally out of scope

## Likely Internal Shape

- sharpen the solver-owned time-state contract in `clode/simulation/_state.py` and the simulator orchestration layer
- add or tighten tests in `test/test_simulation_contracts.py` and any exact-regression coverage that needs to spell out the supported continuation cases
- update `docs/continuation.md` and the `.design` notes so they describe the lived semantics rather than a hoped-for helper
- keep the current OpenCL kernels and single shared-window execution model unless a very small guardrail change proves clearly worthwhile

## Design Constraints

- no public continuation-helper API in this PR
- no device-side per-work-item `t0` redesign in this PR
- no public config redesign in this PR
- preserve the landed solver-state, output-policy, observer-state, observer-parameter, execution-setting, and stepper-definition boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no implicit or IMEX solver implementation in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no broader public continuation-policy API in the same PR
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Make the supported and unsupported continuation cases explicit in the simulator-side state model.
2. Add runtime or contract-test guardrails around shared-window continuation after diverged attained `tf`.
3. Tighten docs and design notes so they describe the honest current semantics.
4. Record the deferred follow-ons explicitly: device-side per-work-item `t0`, broader numerical-helper adoption, adaptive-controller work, and later repo-surface packaging hygiene.

## Code-Facing Checklist

- `clode/simulation/_state.py`: keep solver-owned continuation state honest and inspectable
- `clode/simulation/base.py`: keep continuation helpers aligned with the supported shared-window contract
- `test/test_simulation_contracts.py`: add or tighten continuation-state behavior coverage
- `docs/continuation.md`: describe the representability limit clearly
- `.design/reference/continuation_timebase_note.md`: keep the internal rationale aligned with the live implementation

## Acceptance Criteria

- the package docs and tests make it explicit that exact shared-window continuation requires one common attained `tf`
- the current low-level continuation surface remains stable while the unsupported divergent-time case becomes harder to misread
- the `.design` docs continue to defer public continuation-helper and packaging work until the internal time model is clearer

## Follow-on If This Lands Cleanly

The next high-value follow-ons should be broader numerical-helper adoption into the remaining observers and the stepper time-base path, then deeper numerical time-base work such as structured-time updates. Public continuation helpers, broader config cleanup, and packaging-facing work should still wait until those internal boundaries stop moving.
