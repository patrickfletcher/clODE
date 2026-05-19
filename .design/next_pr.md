# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Single source of truth for execution-setting defaults and compatibility resolution

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- The first-pass solver-state boundary is live: IVP owns next-solve problem data, `clode/simulation/_state.py` owns Python-side solver state and fetched-output caches, and `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners.
- Integration settings and output/storage policy are now explicitly split all the way into the OpenCL layer: common buffers carry integration settings only, trajectory buffers own `max_store`/`nout`, and runtime observer settings no longer include event timestamp capacity.
- The observer-definition and observer-state cleanup is now landed: built-in observer definitions resolve through `ObserverDefinition` and `ResolvedObserverSpec`, and the active kernels plus matched struct lookups now use observer-state naming consistently.
- Public `SolverParams` and `ObserverParams` remain thin compatibility bundles; broader public config redesign is still intentionally deferred.

## Why this should be next

The internal execution path is now cleaner than the compatibility surface above it. The next concrete friction point is execution-setting default resolution.

Today the same settings are still duplicated across `SolverParams`, `_IntegrationSettings`, `_TrajectoryOutputSettings`, and the simulator constructor signatures. Concrete drift already exists in the live code: `SolverParams.dtmax` and `_IntegrationSettings.dtmax` default to `0.5`, while `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` still default `dtmax` to `1.0`; `FeatureSimulator.max_store` still defaults to `10000000`, while `SolverParams.max_store` defaults to `1000000`.

This is the right next PR because it is narrower than a public API redesign, directly addresses predictable config behavior, and gives later stepper-definition work one canonical path for default resolution instead of another moving compatibility layer.

## Scope

- define one canonical internal source for default integration settings and trajectory-output settings
- make `SolverParams` and simulator constructor resolution derive from that source instead of open-coded duplicated defaults
- ensure transient, trajectory, and feature simulators resolve scalar solver arguments and prebuilt `SolverParams` bundles through the same policy
- keep the public constructor signatures and `SolverParams` compatibility bundle intact for now
- document what remains deferred for later stepper-definition cleanup and broader public config redesign

## Likely Internal Shape

- keep `_IntegrationSettings` and `_TrajectoryOutputSettings` as the semantic split between integration policy and output/storage policy
- add explicit canonical defaults or factory helpers in `clode/simulation/params.py`
- make `SolverParams` construction flow from that canonical settings split rather than repeating literal defaults in multiple classes
- let `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` build their initial solver bundle through one shared default-resolution helper
- treat the now-canonical `ObserverParams` defaults as the reference pattern, not as scope to reopen in this PR

## Design Constraints

- no public config API redesign yet
- keep thin compatibility layers close to the current user-facing simulator and parameter bundles until execution-setting and stepper semantics stop moving
- preserve the landed solver-state, output-policy, and observer-state boundaries
- keep fetched outputs and transfer caches as derived data, not semantic owners
- preserve current kernel specialization by precision, stepper, observer, and problem shape unless a narrower path proves clearly better

## Non-goals

- no new public settings classes
- no constructor signature removals or compatibility-breaking renames
- no Python-owned stepper-definition redesign in the same PR
- no public continuation-policy helper in the same PR
- no chunked trajectory streaming or ensemble batching in the same PR
- no multi-device work

## Suggested Implementation Slices

1. Canonicalize default integration and trajectory-output settings in `clode/simulation/params.py`.
2. Route simulator constructor default resolution through one shared helper instead of per-class literals.
3. Add focused regressions proving constructor defaults match the canonical bundles and that scalar-argument resolution matches bundle-based resolution.
4. Record the next follow-on boundary explicitly: stepper-definition cleanup, then broader public config redesign.

## Code-Facing Checklist

- `clode/simulation/params.py`: canonical settings defaults and bundle-conversion helpers
- `clode/simulation/base.py`: shared constructor/default-resolution path
- `clode/simulation/trajectory.py`: trajectory-specific output-policy defaults and invalidation behavior
- `clode/simulation/features.py`: compatibility-safe handling of inherited solver/output arguments without divergent literals
- `test/test_simulation_contracts.py`: constructor-default and default-resolution coverage
- `test/test_opencl_executors.py`: keep runtime invalidation semantics green while settings resolution is refactored

## Acceptance Criteria

- one canonical internal source defines default integration settings and trajectory-output settings
- `SolverParams` defaults and simulator constructor defaults no longer drift
- scalar constructor arguments and prebuilt `SolverParams` bundles resolve to equivalent execution settings
- output-only changes continue to avoid unnecessary runtime invalidation
- the root `.design` docs explicitly state that the next major follow-on is Python-owned stepper-definition cleanup, not immediate public config redesign

## Follow-on If This Lands Cleanly

The next high-value PR should be a cleaner Python-owned stepper-definition model. Public continuation helpers, chunked trajectory and batching work, and any broader public config redesign should follow only once those internal execution-setting and stepper boundaries stop moving.
