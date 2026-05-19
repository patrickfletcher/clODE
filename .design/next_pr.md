# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Broader numerical-helper adoption and time-base groundwork

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, observer-parameter, and stepper-definition boundaries are live.
- `Simulator.advance_tspan_to_attained_final_time()` now provides a narrow exact-continuation helper for the shared-final-time case.
- `FeatureSimulator` now resolves legacy `observer_*` compatibility inputs through one canonical `ObserverParams` path and copies caller-provided bundles instead of aliasing them.
- `clode/kernels/clODE_utilities.cl` now carries a small shared compensated-accumulation helper layer.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- `test/kernel_components/` now provides direct OpenCL component coverage for helper math plus `basic` and `basicall` observer contracts.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The shared helper layer is now real rather than aspirational, but most of the observer and stepper paths still do not use it. That leaves numerical robustness improvements half-landed and the current component-test layer underused.

This is the right next PR because it builds directly on the kernel-helper and continuation groundwork that just landed. It improves single-precision resilience and internal numerical consistency without opening a larger per-work-item time-base redesign too early.

## Scope

- broaden shared compensated or structured helper usage into remaining observer and stepper-time paths where it clearly improves robustness or readability
- add or extend direct component tests where helper adoption changes a contract that is easier to validate below the full simulator stack
- keep the public API stable in this PR
- keep public packaging hygiene, citation metadata, and release-tag cleanup intentionally out of scope

## Likely Internal Shape

- extend `clode/kernels/clODE_utilities.cl` only where the helper layer actually reduces duplicated or fragile kernel logic
- apply those helpers to the remaining observer kernels and any stepper time-base code that still obviously benefits
- extend `test/kernel_components/` and exact-regression coverage only where the helper-adoption pass reveals a contract worth pinning down directly
- leave broader per-work-item `t0` or full structured-time redesign to a later, more explicit PR if it still looks necessary afterward

## Design Constraints

- no device-side per-work-item `t0` redesign in this PR
- no public config redesign in this PR
- preserve the landed solver-state, continuation-helper, output-policy, observer-state, observer-parameter, execution-setting, and stepper-definition boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no implicit or IMEX solver implementation in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no broader continuation-policy redesign in the same PR
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Identify the remaining observer and stepper paths that still duplicate numerically delicate accumulation or time-base logic.
2. Apply the shared helper layer only where it materially clarifies the kernel code or improves robustness.
3. Extend direct component or exact-regression coverage where the changed kernel contract needs tighter pinning.
4. Record the deferred follow-ons explicitly: deeper structured-time work, per-work-item `t0`, adaptive-controller work, and later repo-surface packaging hygiene.

## Code-Facing Checklist

- `clode/kernels/clODE_utilities.cl`: keep the helper layer small, reusable, and numerically explicit
- relevant kernels under `clode/kernels/observers/` and `clode/kernels/steppers/`: adopt the helper layer only where it buys real value
- `test/kernel_components/`: extend direct contracts where helper adoption exposes a better component seam
- `test/core_numerics/`: keep exact-regression coverage honest when helper adoption changes a numerical path
- `.design/reference/testing_audit.md` and `.design/reference/continuation_timebase_note.md`: keep the active rationale aligned with the lived implementation

## Acceptance Criteria

- the remaining observer or stepper paths that obviously benefit now use the shared helper layer instead of duplicated numerically delicate logic
- direct component or exact-regression coverage grows where the helper-adoption pass changes a contract that should stay easy to validate
- the `.design` docs continue to defer deeper structured-time redesign, broader continuation-policy work, and packaging cleanup until the numerical and state-model surface settles further

## Follow-on If This Lands Cleanly

The next high-value follow-ons should be a targeted test-surface audit, then deeper time-base work such as structured-time updates or per-work-item `t0` only if the remaining continuation pressure still justifies it. Broader config cleanup and packaging-facing work should still wait until those internal boundaries stop moving.
