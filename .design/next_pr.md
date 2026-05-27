# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Summary observer family and selected-variable policy

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The observer concept audit is now captured in `.design/reference/observer_concept_audit.md`, including the conclusion that clODE observers are broader than `solve_ivp`-style event functions and that the best first proof target is the summary-observer family.
- The simulation state and output ownership pass is wrapped tightly enough that integration settings, trajectory output policy, observer runtime settings, event-output policy, solver state, persistent observer state, and fetched outputs have clearer named homes on the Python side.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.
- `basic` and `basicall` are one-pass, no-event summary reducers that already share the same online reduction logic and compensated-mean numerics but still exist as separate built-ins with coarse selection scope.
- `basic` tracks one selected feature variable, while `basicall` eagerly tracks all state variables plus all auxiliary variables. There is still no middle ground for "some but not all variables" or for narrower summary groups that would reduce persistent state and output size.
- The heavier event observers remain important, but they should not be the first place where the new declaration model is proven.

## Why this should be next

The audit is complete enough to stop being the active target. The next useful move is to prove the observer declaration model on the simplest family that still matters to users.

This is the right next cut because summary observers already cover an essential clODE workflow: online trajectory reduction without full trajectory storage. They also expose the clearest missing user control today, namely choosing which variables and summary groups to track so state and output size scale with the requested readout rather than with all variables.

It is lower risk than starting with event detectors because it avoids warmup passes, crossing interpolation, and event-retention semantics while still exercising Python-side declaration, layout derivation, build inputs, and compatibility preservation.

## Scope

- unify `basic` and `basicall` behind one internal summary-observer family
- make tracked state variables, tracked auxiliary variables, and slope-summary inclusion explicit parts of the internal declaration model
- preserve current `Observer.basic` and `Observer.basic_all_variables` behavior as compatibility presets over that smaller model
- keep the one-observer-per-build compile-time model, but let summary state layout and feature schema scale with the selected summary scope rather than always one variable or all variables
- keep event detectors and trajectory variable-subset storage out of this PR

## Likely Internal Shape

- introduce one summary-family declaration or adjacent helper layer that can express current one-variable and all-variable summary modes as presets
- derive feature names, persistent layout, and any build- or layout-signature inputs from explicit selected summary groups rather than from a hard-coded observer mode split
- carry selected indices or related summary-family configuration in a narrower internal path than the current broad `ObserverRuntimeSettings` bundle when that is needed to make state size truly depend on selected scope
- keep `_opencl` as the consumer of the resolved summary-family declaration rather than re-specifying summary structure ad hoc

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, solver-owned status boundary, continuation semantics, and current proof layer
- do not reintroduce solver-owned step, status, or time diagnostics through observer outputs or observer-private state
- keep semantic ownership in Python first and treat `_opencl` as the execution consumer rather than the authoritative definition layer
- no user-facing custom-observer DSL, broad inheritance hierarchy, or full code-conversion surface in the same PR
- avoid turning this PR into a benchmark, register-pressure, or "add every missing observer" campaign; the point is a smaller declaration model and better scope control for summary workflows
- preserve current feature names and behavior for the existing public `basic` and `basic_all_variables` presets unless a change is clearly justified and documented
- keep public docs and design notes explicit about what is clarified here and what remains deferred

## Non-goals

- no implicit or IMEX solver work
- no multi-device work
- no broader diverged-time continuation-policy redesign or matched device-side current-time model in the same PR
- no public bundled solver-stats object in the same PR
- no solver-family or stepper-model redesign in the same PR
- no public observer-state fetch API or public rename of `ObserverOutput` in the same PR
- no trajectory variable-subset implementation in the same PR
- no event-observer redesign beyond whatever compatibility hooks are needed to keep the family model coherent
- no citation metadata, release-tag, or broader repo-surface cleanup in the same PR

## Suggested Implementation Slices

1. Introduce the internal summary-family declaration and map current `basic` and `basicall` modes onto it.
2. Make resolved summary feature names and state layout depend on selected state or aux scope and summary groups.
3. Add one narrow selected-subset contract beyond the current one-variable and all-variable presets to prove the model without touching the event observers.

## Code-Facing Checklist

- `clode/observers/_definitions.py`, `clode/observers/types.py`, `clode/observers/metadata.py`: introduce the summary-family declaration and compatibility mappings for current summary modes
- `clode/simulation/features.py`, `clode/simulation/results.py`: keep simulator configuration and `ObserverOutput` coherent with the narrower summary selection model
- `clode/_opencl/observer_metadata.py`, `clode/_opencl/executors.py`, `clode/_opencl/source_builder.py`: make resolved summary scope drive feature metadata, layout shape, and any build-signature changes
- `clode/kernels/observers/observer_basic.clh`, `clode/kernels/observers/observer_basic_allVar.clh`, and any shared replacement helper: keep the summary kernel contract small and explicit
- targeted tests under `test/core_numerics/test_features_basicall.py`, `test/kernel_components/test_kernel_math.py`, `test/test_opencl_models.py`, and `test/test_simulation_contracts.py`: preserve current summary presets and add one narrower selection contract

## Acceptance Criteria

- current `Observer.basic` and `Observer.basic_all_variables` remain supported through compatibility presets over one internal summary family
- at least one narrower selected-summary scope exists beyond the current one-variable or all-variable extremes
- the resolved summary feature schema and persistent layout shrink when the selected summary scope shrinks
- the updated code and docs make the later extension of the same selection model to heavier observers more explicit than it is today

## Follow-on If This Lands Cleanly

1. decide which parts of the same selection model should extend to local-extrema, threshold, or neighborhood observers
2. add a small set of additional dynamical-systems-oriented built-in observers or features on top of the smaller declaration model
3. evaluate whether any user-authored or code-converted observer surface is justified once the built-in declaration model stops moving
