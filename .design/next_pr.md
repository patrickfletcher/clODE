# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

InitialValueProblem first, with built-in batch semantics at the simulator boundary

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- Split-window continuation regressions now pass on the live PyOpenCL path, with solver-owned absolute time and persisted stochastic continuation state.
- The public simulator API already requires default parameter values and default initial-state values through the `parameters` and `variables` mappings.
- Root flat modules such as `clode.solver` and `clode.features` are compatibility barrels only; new work should target canonical packages unless the task is explicitly about compatibility cleanup.
- Archived migration notes remain useful for reproduction details, but they are not the source of truth for the live package layout.

## Why this should be next

Continuation correctness is no longer the blocker. The live split-window regressions in `test/core_numerics/test_features_basicall.py` and `test/core_numerics/test_stochastic.py` now pass on the stable workspace runtime (`CLODE_TEST_PLATFORM_ID=0`, `CLODE_TEST_DEVICE_ID=0`).

The next gap is at the user-facing problem and batch-input boundary:

- the public API already requires default parameter values and default initial-state values, but those defaults still live as raw `Simulator` constructor mappings rather than one explicit `InitialValueProblem`
- `set_ensemble()` and `set_repeat_ensemble()` already expose a distinct ensemble concern, but that concern is still encoded as simulator-side array shaping and broadcasting helpers
- simulator classes therefore still mix IVP/default construction, ensemble shaping, and solve orchestration in one place

The likely first move is to let an IVP cover both the usual size-`(1,)` case and basic batched parameter or initial-state inputs, with shape metadata and helper functions where useful. A separate public `IVPEnsemble` class should be deferred until the live code shows behavior that is meaningfully richer than normalized arrays plus remembered shape.

## Scope

- define simulator classes as orchestration objects that compose an `InitialValueProblem`, the `_opencl` executor/runtime, and observer or trajectory policy
- introduce an `InitialValueProblem` semantic model centered on RHS plus default parameter values and default initial-state values
- let that IVP model own the size-`(1,)` case and basic batched initial-state or parameter semantics, including remembered ensemble shape where that improves result reshaping or plotting workflows
- add helper functions for grid, random, quasi-random, and repeat-style batch generation around the IVP model before committing to a separate public ensemble class
- keep lower-level shared definition metadata such as `ProblemInfo`, `ProblemShape`, and `RhsSource` as derived or internal support concepts unless a separate user-facing layer proves necessary later
- keep lower-level solver-state, observer-state, and stepper-definition cleanup as follow-on work unless a small internal adapter falls out naturally

## Non-goals

- no explicit per-work-item solver-state redesign in this pass
- no observer-definition or `ObserverData` redesign
- no stepper-definition model work
- no chunked trajectory streaming or storage redesign
- no implicit solver addition
- no multi-device work
- no Random123 adoption yet
- no broad public API redesign
- no broad package, module, or kernel-asset relocation

## Acceptance Criteria

- split-window continuation regressions remain green for `basicall` features and seeded stochastic Euler
- simulator constructors can be understood as building or consuming an explicit `InitialValueProblem` rather than directly owning the default-state semantics themselves
- simulator classes are more obviously orchestration-focused, with a clearer internal home for batch shaping and remembered result shape instead of keeping ensemble generation primarily in `Simulator` helper methods
- no separate public `ProblemDefinition` layer is required for the first semantic pass; shared static definition metadata can remain derived or internal
- a separate public `IVPEnsemble` class is not required in the first semantic pass; batching can remain IVP-owned or helper-driven unless richer behavior proves necessary
- lower-level solver-state, observer-state, and stepper-definition cleanup is not required to land this pass beyond incidental adapter changes

## Follow-on If This Lands Cleanly

The next high-value steps should be an explicit per-work-item solver-state model and public continuation-policy helpers, then separation of integration state from output and storage policy. A dedicated ensemble class can be revisited later if IVP-owned batching and helper functions prove too limiting.
