# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Adaptive-time compensated time base and observer elapsed-time follow-on

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state, observer-definition, observer-parameter, and stepper-definition boundaries are live.
- `Simulator.advance_tspan_to_attained_final_time()` now provides a narrow exact-continuation helper for the shared-final-time case.
- `FeatureSimulator` now resolves legacy `observer_*` compatibility inputs through one canonical `ObserverParams` path and copies caller-provided bundles instead of aliasing them.
- `clode/kernels/clODE_utilities.cl` now carries the shared helper layer for compensated means, compensated time pairs, fixed-step counter time reconstruction, threshold timestamps, and bounded three-sample extrema.
- The `basic` and `basicall` observers now use compensated integral accumulation for their time-weighted means.
- Fixed-step live steppers now reconstruct absolute time from `t0 + step * dt` with a 64-bit step counter.
- `threshold_2` now uses inverse-linear timestamps for stored up/down threshold transitions, and `local_max` now uses the shared bounded three-sample max/min helpers.
- `test/kernel_components/` now provides direct OpenCL component coverage for the helper prototypes plus the `basic`, `basicall`, `threshold_2`, and `local_max` observer contracts touched by the numerics work.
- Adaptive steppers still derive the next absolute time from one float32 absolute-time value plus `dt`, and several observer elapsed-time or duration paths still depend on subtracting large float32 absolutes.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The conservative live-helper slice is now landed. The largest remaining numerics risk with clear evidence is the adaptive-time path: long absolute-time runs and large-origin feature windows still depend on one float32 absolute-time value and on subtracting large float32 timestamps in observer bookkeeping.

This keeps the package aligned with the main product story: fast large-ensemble workflows where users mostly want final states or on-device features, and where single-precision robustness needs to be demonstrably better rather than assumed.

This is the right next PR because it attacks the remaining root cause instead of broadening helper adoption opportunistically. A solver-owned compensated adaptive-time path plus observer-facing relative elapsed bookkeeping should improve real behavior where the fixed-step counter shortcut does not apply, while keeping the public API stable and leaving implicit-method work for later.

## Scope

- implement a solver-owned dual-realtype compensated time path for adaptive steppers and any shared adaptive-time bookkeeping that still relies on one float32 absolute time plus `dt`
- provide observer-facing relative elapsed-time bookkeeping, or an equivalent derived channel, for the time-weighted and duration-sensitive paths that currently subtract large float32 absolutes
- thread the richer adaptive-time bookkeeping through the relevant kernel structs, stepper interfaces, and observer call sites without widening the public compatibility surface
- extend component and focused simulator coverage for large-origin adaptive time and elapsed-time cases
- update docs and design notes so they distinguish the landed fixed-step or observer helper work from the still-open adaptive-time follow-on
- keep the public API stable in this PR
- keep broader helper rollout, public packaging hygiene, citation metadata, and release-tag cleanup intentionally out of scope

## Likely Internal Shape

- use the existing compensated-time helper primitives as the implementation home for the adaptive solver-owned time pair rather than duplicating time logic across steppers
- make adaptive steppers consume a richer absolute-time representation while keeping their numerical policy explicit and inspectable
- derive observer elapsed-time and duration bookkeeping from the richer time base, or from a solver-owned relative elapsed channel, where large-origin subtraction is the real failure mode
- update the component tests so they cover the adopted adaptive-time behavior in addition to the already-landed helper and observer contracts
- keep the public docs explicit about which time paths are now live and which broader observer rollouts remain deferred

## Design Constraints

- no public config redesign in this PR
- preserve the landed solver-state, continuation-helper, output-policy, observer-state, observer-parameter, execution-setting, and stepper-definition boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- do not disturb the landed fixed-step counter-reconstructed time path while adding the adaptive-time follow-on
- do not treat the adaptive-time follow-on as justification for broad observer-helper rollout; any additional observer changes should be limited to elapsed-time bookkeeping that the new time representation makes necessary
- no helper in this PR should depend on a general nonlinear system solver; anything that needs that should still be deferred with the later implicit-method work
- public docs should continue recommending autonomous feature windows near `t = 0` when that does not change model semantics
- public docs should distinguish the landed fixed-step and observer helper adoption from the remaining adaptive-time work instead of implying everything is already solved
- public docs should compare current package features to relevant algorithmic alternatives, not narrate the package's development history
- leave citation metadata, release-tag alignment, and similar repo-surface packaging hygiene for later

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no broad helper rollout across all remaining observers in the same PR
- no implicit or IMEX solver implementation in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work
- no broader diverged-time continuation-policy redesign in the same PR
- no citation metadata or release-tag cleanup in the same PR

## Suggested Implementation Slices

1. Introduce the adaptive solver-owned compensated time representation and wire it through the adaptive stepper entrypoints.
2. Move the observer elapsed-time and duration-sensitive paths that still subtract large float32 absolutes onto the richer time bookkeeping.
3. Extend component and focused simulator tests for large-origin adaptive-time and observer-elapsed-time cases.
4. Update the docs and design notes to reflect the new adaptive-time behavior and the remaining broader deferrals.

## Code-Facing Checklist

- `clode/kernels/clODE_utilities.cl`: keep the shared helper layer as the only implementation home for compensated adaptive-time primitives
- `clode/kernels/steppers/` plus the adaptive driver entrypoints: move adaptive absolute time off one-float absolute-time bookkeeping
- `clode/kernels/observers/*.clh`: limit changes to elapsed-time or duration bookkeeping that depends on subtracting large float32 absolutes
- `clode/kernels/clODE_struct_defs.cl` and the host-side matched structs: thread any richer adaptive-time state through the kernel boundary without reopening the public compatibility layer
- `test/kernel_components/test_kernel_math.py`: keep direct helper coverage alongside the already-landed live observer checks and any new adaptive-time component contracts
- `docs/numerical_accuracy.md`: explain which time-base fixes are live for fixed-step, which are live for adaptive paths after this PR, and which broader work remains deferred
- `.design/reference/single_precision_numerics_note.md` and `.design/reference/continuation_timebase_note.md`: keep the active rationale aligned with the lived implementation

## Acceptance Criteria

- the adaptive live solver path no longer relies solely on one float32 absolute-time value plus `dt`
- the observer paths that are updated for elapsed-time or duration bookkeeping no longer rely solely on subtracting large float32 absolutes where the richer time representation is available
- component tests and focused simulator tests cover the adopted adaptive-time behavior in addition to the retained helper and observer contracts
- the public docs and design notes clearly distinguish the landed fixed-step or observer helper behavior from the new adaptive-time changes and from the still-deferred broader rollout or implicit-method work
- the `.design` docs continue deferring broader helper rollout, diverged-time continuation redesign, implicit-solver or nonlinear-solver work, and packaging cleanup until the numerical and state-model surface settles further

## Follow-on If This Lands Cleanly

If this lands cleanly, the next grouped follow-ons should stay selective and ordered:

1. execution-state and invalidation hardening: solver-state or solution-buffer abstraction, explicit rebuild-policy boundaries, and direct build-key or component coverage
2. numerical validation and evidence: a slim exact-solution suite plus public demonstrations of where the adopted precision safeguards matter
3. large-ensemble ergonomics: IVP-side batch-generation helpers and device-capacity ensemble batching before broader trajectory-output expansion

Keep broader observer-helper rollout, source-assembly reshaping, broader PyOpenCL helper leverage, and any helper that would need general nonlinear system solve machinery deferred until those grouped follow-ons settle.
