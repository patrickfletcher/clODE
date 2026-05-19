# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Explicit observer-definition model and observer runtime/state cleanup

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state boundary is live: IVP owns next-solve problem data, `clode/simulation/_state.py` owns Python-side solver state and fetched-output caches, and `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners.
- Integration settings and output/storage policy are now explicitly split all the way into the OpenCL layer: common buffers carry integration settings only, trajectory buffers own `max_store`/`nout`, and runtime observer settings no longer include event timestamp capacity.
- Public `SolverParams` and `ObserverParams` remain thin compatibility bundles; broader public config redesign is still intentionally deferred.

## Why this should be next

The solver-state cleanup and output/storage split both landed cleanly. The next structural friction point is now concentrated in the observer model.

Today feature-name generation, warmup behavior, persistent observer state layout, optional event-output layout, and runtime observer settings are still spread across `clode/observers/metadata.py`, `clode/_opencl/observer_metadata.py`, `_opencl/executors.py`, and the observer kernel include tree. That makes observer work harder to reason about than the now-stable solver/output boundary beneath it.

This is the right next PR because it builds directly on the newly landed boundaries, keeps public compatibility layers thin and close to the existing user-facing API, and establishes a much cleaner prerequisite for any later public config redesign.

## Scope

- introduce one explicit internal observer-definition layer that owns feature names, warmup/two-pass requirements, persistent state layout, runtime settings, and optional event-output layout boundaries
- keep persistent observer state distinct from solver state and from optional event-output policy
- make layout-affecting observer changes explicit so rebuild/reallocation rules follow clear policy boundaries
- preserve the current public `FeatureSimulator` and `ObserverParams` surface except for incidental compatibility-safe polish
- document which internal cleanups should finish before any broader public config redesign begins

## Likely Internal Shape

- keep the public `Observer` enum and `ObserverParams` compatibility bundle as the user-facing selector and parameter surface for now
- add one Python-owned `ObserverDefinition` per built-in observer that owns stable semantic facts: observer id/name, warmup/two-pass requirement, feature-schema factory, runtime-setting interpretation, and the distinction between persistent state and optional event-output layout
- resolve a selected observer plus `ProblemInfo`, runtime settings, precision, and event-output settings into a problem-shaped `ResolvedObserverSpec` or equivalent value object used by `_opencl`
- let that resolved spec provide the inputs currently scattered across `clode/observers/metadata.py`, `clode/_opencl/observer_metadata.py`, `_opencl/executors.py`, and `_opencl/registry.py`
- keep kernel implementations in the existing observer `.clh` files for now, but make build defines and layout decisions flow from the Python observer definition instead of parallel hardcoded conditionals
- keep `ObserverRuntimeSettings` and `EventOutputSettings` as separate concerns: runtime tuning versus layout/output policy

## Design Constraints

- no public observer or config API redesign yet
- keep thin compatibility layers close to the current user-facing simulator and parameter bundles until observer and stepper semantics stop moving
- keep event timestamp capacity as explicit output/layout policy rather than generic runtime observer config
- preserve the landed solver-state and output-policy boundary; do not move observer state into solver state
- keep fetched outputs and transfer caches as derived data, not semantic owners
- preserve current kernel specialization by precision, observer, and problem shape unless a narrower path proves clearly better

## Non-goals

- no public custom-observer API
- no observer-specific public parameter classes yet
- no chunked trajectory streaming or ensemble batching in the same PR
- no Python-owned stepper-definition redesign in the same PR
- no public continuation-policy helper in the same PR
- no multi-device work

## Suggested Implementation Slices

1. Introduce explicit Python-owned observer definitions and route feature-name generation plus warmup/two-pass facts through them.
2. Separate runtime observer settings, persistent observer state layout, and optional event-output layout in the metadata and executor paths.
3. Narrow program rebuild and buffer reallocation rules so only explicit layout changes force them.
4. Add focused regressions for runtime-only observer updates, event-layout changes, and continued numerical correctness.
5. Record the prerequisites for later public config redesign once the internal observer and stepper boundaries are stable enough.

## Code-Facing Checklist

- `clode/observers/metadata.py`: centralize feature-name and observer-definition facts behind a clearer internal model.
- `clode/observers/types.py`: keep the public compatibility bundle intact while exposing only the internal split views needed by the runtime.
- `clode/_opencl/observer_metadata.py`: separate persistent observer-state layout from optional event-output layout more explicitly.
- `clode/_opencl/executors.py`: make observer rebuild/reallocation logic follow explicit runtime-vs-layout distinctions.
- `clode/simulation/features.py`: preserve the current user-facing surface while keeping cache invalidation aligned with the refined observer boundaries.
- `clode/kernels/observers.cl` and related observer include files: touch only where needed to align the runtime/layout split and warmup semantics.
- `test/test_opencl_executors.py` and `test/test_simulation_contracts.py`: add or adjust coverage for runtime-only observer updates versus layout-changing updates.
- `test/core_numerics/test_features_basicall.py` and `test/core_numerics/test_stochastic.py`: keep continuation and observer correctness green while the cleanup lands.

## Acceptance Criteria

- there is one explicit internal home for observer definitions, runtime settings, persistent observer-state layout, optional event-output layout, feature names, and warmup requirements
- runtime-only observer updates no longer trigger unnecessary rebuilds or buffer reallocations
- optional event-output capacity remains explicit layout policy rather than generic runtime observer configuration
- the landed solver-state and output-policy boundary remains intact
- the current public `FeatureSimulator` and `ObserverParams` compatibility paths remain green unless a deliberate, documented change is made
- the root `.design` docs explicitly state that broader public config redesign should wait until observer-definition and stepper-definition cleanup land

## Follow-on If This Lands Cleanly

The next few high-value PRs should be a cleaner Python-owned stepper-definition model, then chunked trajectory/batching helpers built on the stabilized solver/output/observer boundaries, and only then a broader public config redesign. Public continuation helpers and richer IVP batch helpers should follow once those internal boundaries stop moving.
