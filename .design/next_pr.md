# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Kernel-math helper foundation and component-test spine

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state boundary is live: IVP owns next-solve problem data, `clode/simulation/_state.py` owns Python-side solver state and fetched-output caches, and `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners.
- Integration settings and output/storage policy are now explicitly split all the way into the OpenCL layer, and execution-setting defaults now resolve through one canonical solver-settings helper.
- The observer-definition and observer-state cleanup is now landed: built-in observer definitions resolve through `ObserverDefinition` and `ResolvedObserverSpec`, and the active kernels plus matched struct lookups now use observer-state naming consistently.
- Built-in stepper traits and OpenCL build mapping now resolve through a Python-owned stepper-definition catalog.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The current codebase is strongest when its internals are explicit and testable. The kernels already contain TODOs around Kahan or TwoSum-style accumulation, FMA choices, and time-accumulation precision, and the live tests still jump from full end-to-end numerics to smaller host-side `_opencl` support checks without a dedicated middle layer.

This is the right next PR because it improves soundness, internal reasoning, and future extensibility for observers and steppers without forcing more public API churn. It also lines up with the current user priority: tighten internals first, then return to broader ergonomics once the numerical and testing substrate is stronger.

## Scope

- introduce or tighten one small shared OpenCL numerical-helper layer for reusable precision-sensitive utilities used by observers and later steppers
- add tiny synthetic component tests around those helpers, build specialization, and observer or stepper contracts that are currently hard to validate without end-to-end runs
- keep kernel specialization, `KernelKind`, and source assembly explicit and legible; audit those boundaries only where the helper and test work touches them
- keep the public API unchanged
- document the deferred follow-ons clearly: continuation-state guardrails, observer-surface cleanup, and deeper numerical time-base work

## Likely Internal Shape

- keep shared numerical helpers in `clode/kernels/clODE_utilities.cl` and/or `clode/kernels/realtype.cl` rather than scattering them ad hoc across observer and stepper kernels
- adopt those helpers only where they actually simplify duplicated or fragile arithmetic in the current kernels
- add a small component-test layer using tiny synthetic models or kernels so observer accumulation and stepper-helper behavior can be checked directly
- keep kernel files separate from Python definition catalogs; do not move OpenCL kernels into Python string literals in this PR

## Design Constraints

- no public API redesign in this PR
- keep thin compatibility layers close to the current user-facing simulator and parameter bundles until internal observer, stepper, and continuation semantics stop moving
- preserve the landed solver-state, output-policy, observer-state, execution-setting, and stepper-definition boundaries
- keep kernel specialization explicit and inspectable rather than hiding it behind a larger meta-build framework
- do not bulk-relocate kernel files or absorb them into Python definitions without strong evidence that it improves reasoning instead of just moving complexity around
- keep fetched outputs and transfer caches as derived data, not semantic owners
- preserve current kernel specialization by precision, stepper, observer, and problem shape unless a narrower path proves clearly better

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no device-side per-work-item `t0` redesign in the same PR
- no implicit or IMEX solver implementation in the same PR
- no continuation helper or broader continuation-policy API in the same PR
- no public config redesign in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work

## Suggested Implementation Slices

1. Audit the current kernels for precision-sensitive duplicated arithmetic and choose the smallest helper set that materially improves clarity or robustness.
2. Add focused component tests for one observer path and one stepper or numerical-helper path without relying only on end-to-end simulator runs.
3. Adopt the helper layer in a small, high-value slice rather than rewriting every observer or stepper kernel at once.
4. Record the deferred follow-ons explicitly: continuation-state guardrails, observer-surface cleanup, adaptive-controller work, and deeper time-base refinements.

## Code-Facing Checklist

- `clode/kernels/clODE_utilities.cl` and `clode/kernels/realtype.cl`: establish the reusable helper boundary clearly
- relevant observer and stepper kernels under `clode/kernels/observers/` and `clode/kernels/steppers/`: adopt helpers only where they buy clarity or robustness
- `test/`: add a component-test layer or equivalent focused tests for the new helper and kernel contracts
- `.design/reference/testing_audit.md`: keep the test-strategy note aligned with the new middle layer
- `.design/reference/pyopencl_leverage_audit.md`: keep the PyOpenCL boundary explicit if helper or test work changes how much host-side infrastructure clODE owns

## Acceptance Criteria

- one small shared numerical-helper layer exists or the current helper boundary is explicitly tightened instead of remaining scattered TODOs
- the test surface now includes a focused middle layer between end-to-end numerics and small host-side support tests
- at least one observer path and one stepper or numerical-helper path become easier to reason about through direct tests or clearer shared utilities
- the root `.design` docs explicitly state that broader public config redesign still waits until these internal numerical and testing boundaries stabilize

## Follow-on If This Lands Cleanly

The next high-value follow-ons should be observer-system and stepper-extension cleanup, continuation-state guardrails, then deeper numerical time-base work such as structured-time updates. Broader public continuation or config API work should still wait until those internal boundaries stop moving.
