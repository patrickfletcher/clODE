# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer-system audit and compatibility-surface cleanup

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, and stepper-definition boundaries are live.
- `clode/kernels/clODE_utilities.cl` now carries a small shared compensated-accumulation helper layer.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- `test/kernel_components/` now provides direct OpenCL component coverage for helper math and a basic-observer contract.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The observer-definition catalog is now good enough that the remaining observer complexity is easier to see clearly. The current package still carries a large legacy `observer_*` scalar surface in `FeatureSimulator`, and some of the built-in observer helper or kernel structure still reads as historical accumulation rather than deliberate design.

This is the right next PR because it targets one of clODE's most distinctive workflows directly: on-device feature extraction. It improves internal reasoning and user-facing ergonomics without jumping into public packaging work or broader API redesign too early.

## Scope

- audit the built-in observer helpers and kernels against the landed `ObserverDefinition` model and remove or simplify complexity that no longer carries its weight
- keep `ObserverParams` as the semantic owner of observer configuration while reducing the amount of observer meaning carried by constructor-level scalar compatibility plumbing
- extend component tests where needed so observer behavior is easier to validate directly than it is today
- keep the public API stable in this PR
- document that citation metadata, release-tag hygiene, and similar public packaging work remain intentionally later than core package-quality work

## Likely Internal Shape

- keep `ObserverDefinition` and `ResolvedObserverSpec` as the semantic boundary for built-ins instead of reopening observer meaning in `_opencl`
- audit `FeatureSimulator` parameter resolution so scalar `observer_*` arguments behave as thin compatibility inputs layered over `ObserverParams`
- extend `test/kernel_components/` only where the observer audit reveals contracts that are still too implicit
- keep kernel files separate from Python definition catalogs; do not move OpenCL kernels into Python string literals in this PR

## Design Constraints

- no public API redesign in this PR
- keep thin compatibility layers close to the current user-facing simulator and parameter bundles until internal observer, stepper, and continuation semantics stop moving
- preserve the landed solver-state, output-policy, observer-state, execution-setting, and stepper-definition boundaries
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- do not bulk-relocate kernel files or absorb them into Python definitions without strong evidence that it improves reasoning instead of just moving complexity around
- keep fetched outputs and transfer caches as derived data, not semantic owners
- preserve current kernel specialization by precision, stepper, observer, and problem shape unless a narrower path proves clearly better
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no device-side per-work-item `t0` redesign in the same PR
- no implicit or IMEX solver implementation in the same PR
- no continuation helper or broader continuation-policy API in the same PR
- no public config redesign in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Audit the built-in observer catalog and kernels for duplication, stale complexity, or unclear semantics.
2. Tighten `FeatureSimulator` observer-parameter resolution so `ObserverParams` is clearly the semantic owner and scalar observer args remain only compatibility shims.
3. Add or extend component tests for any audited observer contract that is still hard to reason about directly.
4. Record the deferred follow-ons explicitly: continuation-state guardrails, broader helper adoption, adaptive-controller work, and later repo-surface packaging hygiene.

## Code-Facing Checklist

- `clode/simulation/features.py`: keep observer parameter resolution and compatibility inputs honest and narrow
- `clode/observers/_definitions.py` and related metadata helpers: keep the semantic observer catalog authoritative
- relevant observer kernels under `clode/kernels/observers/`: simplify only where the audit shows clear value
- `test/kernel_components/`: extend direct observer contract coverage where it buys clarity
- `.design/reference/testing_audit.md`: keep the component-test note aligned with the lived test surface

## Acceptance Criteria

- the observer catalog and built-in observer kernels are easier to reason about after one explicit audit or cleanup pass
- `FeatureSimulator` clearly treats `ObserverParams` as the semantic owner and scalar observer arguments as compatibility inputs
- the component-test layer covers at least one additional observer contract if the audit exposes a still-implicit behavior boundary
- the root `.design` docs explicitly state that public packaging hygiene remains later than the internal package-quality work

## Follow-on If This Lands Cleanly

The next high-value follow-ons should be continuation-state guardrails, broader numerical-helper adoption into remaining observers and the stepper time-base path, then deeper numerical time-base work such as structured-time updates. Broader public continuation or config API work should still wait until those internal boundaries stop moving.
