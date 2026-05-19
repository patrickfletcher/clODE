# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Python-owned stepper-definition model

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state boundary is live: IVP owns next-solve problem data, `clode/simulation/_state.py` owns Python-side solver state and fetched-output caches, and `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners.
- Integration settings and output/storage policy are now explicitly split all the way into the OpenCL layer, and execution-setting defaults now resolve through one canonical solver-settings helper.
- The observer-definition and observer-state cleanup is now landed: built-in observer definitions resolve through `ObserverDefinition` and `ResolvedObserverSpec`, and the active kernels plus matched struct lookups now use observer-state naming consistently.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still intentionally deferred.

## Why this should be next

The remaining structural friction is stepper semantics. Right now method traits and build mapping are still split across the public `Stepper` enum, `_opencl/registry.py`, source building, and the kernel include tree.

This is the right next PR because the execution-setting boundary is now stable enough that stepper work no longer has to clean up drifting defaults at the same time. A Python-owned stepper-definition layer should make source assembly, later implicit-stepper work, and eventual public config redesign much easier to reason about.

## Scope

- introduce one internal stepper-definition layer that owns semantic traits and OpenCL build mapping for the built-in steppers
- keep the public `Stepper` enum as the compatibility selector for now
- route source-building and runtime validation through that internal definition layer instead of raw string tables
- make fixed versus adaptive and deterministic versus stochastic traits explicit in Python rather than implicit in scattered kernel or registry logic
- document what remains deferred for continuation helpers, implicit steppers, and broader public config redesign

## Likely Internal Shape

- keep the public `Stepper` enum as the selector used by simulators and callers for now
- add a `StepperDefinition` or equivalent internal value object that owns stable semantic facts such as method name, build define, fixed versus adaptive, deterministic versus stochastic, and future explicit versus implicit traits
- resolve the selected stepper into one internal object consumed by `_opencl/registry.py`, `_opencl/source_builder.py`, and executor construction
- keep the current kernel implementations in `clode/kernels/steppers.cl` and `clode/kernels/steppers/*.clh` for now, but make selection and mapping flow from Python definitions rather than parallel string conditionals

## Design Constraints

- no public config API redesign yet
- keep thin compatibility layers close to the current user-facing simulator and parameter bundles until stepper semantics stop moving
- preserve the landed solver-state, output-policy, observer-state, and execution-setting boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- preserve current kernel specialization by precision, stepper, observer, and problem shape unless a narrower path proves clearly better

## Non-goals

- no constructor signature removals or compatibility-breaking renames
- no implicit or IMEX solver implementation in the same PR
- no public continuation-policy helper in the same PR
- no public config redesign in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work

## Suggested Implementation Slices

1. Introduce explicit Python-owned stepper definitions and route validation and build-define lookup through them.
2. Make fixed/adaptive and deterministic/stochastic traits explicit in the internal model.
3. Update source building and executor construction to consume the new stepper-definition layer.
4. Add focused regressions for stepper-definition lookup, source-builder options, and compatibility with the current public enum surface.
5. Record the follow-on boundary explicitly: continuation helpers, implicit-stepper groundwork, and only later broader public config redesign.

## Code-Facing Checklist

- `clode/simulation/base.py`: keep the public `Stepper` compatibility path intact while narrowing its semantic job
- `clode/_opencl/registry.py`: replace raw stepper string tables with stepper-definition-driven lookup
- `clode/_opencl/source_builder.py`: consume the new stepper-definition layer when producing build options
- `clode/_opencl/models.py`: keep build-key semantics aligned with the refined stepper mapping
- `clode/kernels/steppers.cl`: touch only where needed to align build defines and method families with the new Python-owned definitions
- `test/test_opencl_models.py`, `test/test_opencl_source_builder.py`, and `test/test_opencl_executors.py`: add or adjust coverage for stepper-definition lookup and source-assembly behavior

## Acceptance Criteria

- one internal stepper-definition layer owns the built-in stepper traits and OpenCL build mapping
- source builder and runtime validation no longer depend on scattered raw stepper string tables
- the current public `Stepper` enum path remains green unless a deliberate, documented change is made
- the root `.design` docs explicitly state that broader public config redesign still waits until stepper semantics stabilize

## Follow-on If This Lands Cleanly

The next high-value follow-ons should be public continuation helpers and implicit-stepper groundwork, then chunked trajectory and batching work, and only then a broader public config redesign once those internal execution and stepper boundaries stop moving.
