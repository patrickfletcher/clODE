# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

InitialValueProblem first, with built-in batch semantics at the simulator boundary

## Current Status

The core semantic pass is now on the branch: `InitialValueProblem` owns default state, default parameters, batch shaping, and Python-backed SciPy-style callability, while simulators can now consume an explicit `ivp=` argument and continue to support the compatibility constructor path.

The remaining work in this PR should stay narrow:

- public-surface cleanup so `clode.problem` and the curated docs center `InitialValueProblem`
- low-risk semantics cleanup where lower-level helper types such as `ProblemInfo` and `RhsSource` are treated as derived support concepts rather than promoted user-facing API
- docs and test polishing tied directly to the landed IVP behavior

Broader batch-generation helpers and generic cross-source SciPy interop should remain follow-on work.

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
- stage the implementation in three passes: settle the IVP API and test matrix first, land IVP ownership and simulator delegation second, and add Python-backed callability third
- introduce an `InitialValueProblem` semantic model centered on RHS plus default parameter values and default initial-state values
- let that IVP model own the size-`(1,)` case and basic batched initial-state or parameter semantics, including remembered ensemble shape where that improves result reshaping or plotting workflows
- when an IVP is built from a Python RHS, preserve that Python callable on the IVP and evaluate a small adapter or `__call__` path that can satisfy SciPy `solve_ivp` without changing the OpenCL `getRHS` contract
- add helper functions for grid, random, quasi-random, and repeat-style batch generation around the IVP model before committing to a separate public ensemble class
- keep lower-level shared definition metadata such as `ProblemInfo`, `ProblemShape`, and `RhsSource` as derived or internal support concepts unless a separate user-facing layer proves necessary later
- keep lower-level solver-state, observer-state, and stepper-definition cleanup as follow-on work unless a small internal adapter falls out naturally

## Non-goals

- no explicit per-work-item solver-state redesign in this pass
- no observer-definition or `ObserverData` redesign
- no stepper-definition model work
- no chunked trajectory streaming or storage redesign
- no implicit solver addition
- no requirement that every IVP source be callable from SciPy in this pass; generic interop for OpenCL-only or XPP-defined problems should remain deferred until the project has a stronger RHS representation than today's converter output and source text
- no multi-device work
- no Random123 adoption yet
- no broad public API redesign
- no broad package, module, or kernel-asset relocation

## Acceptance Criteria

- split-window continuation regressions remain green for `basicall` features and seeded stochastic Euler
- the API and test plan in `.design/reference/ivp_api_test_plan.md` remains consistent with the landed implementation, or is updated in the same PR when implementation evidence changes the plan
- simulator constructors can be understood as building or consuming an explicit `InitialValueProblem` rather than directly owning the default-state semantics themselves
- simulator classes are more obviously orchestration-focused, with a clearer internal home for batch shaping and remembered result shape instead of keeping ensemble generation primarily in `Simulator` helper methods
- if SciPy callability lands in this pass, it does so as a Python-backed IVP adapter rather than by reverse-mapping OpenCL source or broadening the converter into a bidirectional IR project
- SciPy is used only as an optional development and test dependency in this pass, and the new tests stay slim and contract-focused rather than trying to freeze the whole evolving API surface
- no separate public `ProblemDefinition` layer is required for the first semantic pass; shared static definition metadata can remain derived or internal
- a separate public `IVPEnsemble` class is not required in the first semantic pass; batching can remain IVP-owned or helper-driven unless richer behavior proves necessary
- lower-level solver-state, observer-state, and stepper-definition cleanup is not required to land this pass beyond incidental adapter changes

## Follow-on If This Lands Cleanly

The next high-value steps should be an explicit per-work-item solver-state model and public continuation-policy helpers, then separation of integration state from output and storage policy. A dedicated ensemble class can be revisited later if IVP-owned batching and helper functions prove too limiting. Broader SciPy or cross-solver interop for OpenCL-only and XPP-defined problems should remain follow-on work coupled to a stronger internal RHS representation rather than to the first IVP semantic pass.
