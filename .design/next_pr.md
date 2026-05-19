# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Explicit solver-state ownership and cache cleanup at the simulator/runtime boundary

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- Split-window continuation regressions now pass on the live PyOpenCL path, with solver-owned absolute time and persisted stochastic continuation state.
- The public simulator API already requires default parameter values and default initial-state values through the `parameters` and `variables` mappings.
- Root flat modules such as `clode.solver` and `clode.features` are compatibility barrels only; new work should target canonical packages unless the task is explicitly about compatibility cleanup.
- Archived migration notes remain useful for reproduction details, but they are not the source of truth for the live package layout.

## Why this should be next

The IVP semantic pass is now in place. `InitialValueProblem` owns default state, default parameters, batch shaping, and Python-backed SciPy-style callability, and the targeted IVP plus regression slices are green.

That work exposed the next real friction point more clearly: ownership of solver-related state and cached host-side data.

Today the runtime semantics are spread across several layers:

- IVP-owned problem data for the next solve
- persistent observer state and observer-layout buffers for feature runs
- simulator-side cached arrays such as `_device_initial_state`, `_device_parameters`, `_device_final_state`, `_device_dt`, and `_device_tf`
- executor-side host mirrors such as `_x0_host`, `_pars_host`, `_xf_host`, `_dt_host`, and `_tf_host`
- device buffers that remain the real execution source of truth during solves

This now feels like the highest-value next cleanup because it creates unnecessary synchronization logic, blurs single-source-of-truth boundaries, and makes continuation and invalidation rules harder to reason about than they need to be.

It is also the cleanup most likely to influence the next observer-definition and output-policy work. If the solver-state boundary stays fuzzy, later feature or trajectory refactors will either duplicate more host/device state or accidentally bake today’s observer and storage assumptions deeper into the runtime contract.

## Scope

- introduce an explicit internal solver-state model with a clear home for requested time window, per-work-item current-start time `t0`, attained final time, current `dt`, and continuation-specific RNG state, with only the minimal status or sync markers needed for this first pass
- define a sharper ownership split between IVP-owned problem data, solver-owned execution state, persistent observer state, and fetched output caches
- reduce duplicated host-side mirrors and synchronization paths between `Simulator` and `_opencl/executors.py`
- make invalidation rules explicit when callers change IVP data, solver parameters, `t_span`, or observer/runtime configuration
- separate compile-time build-spec decisions from runtime solver or observer state so rebuilds follow explicit kernel-specialization inputs rather than incidental cache invalidation
- keep later time-chunking and ensemble-batching as orchestration policies layered on explicit solver state and output policy rather than as hidden side effects of `max_steps`, `max_store`, or simulator caches
- leave observer metadata and feature-buffer code on the current public model, but make the solver-state boundary explicit enough that later observer-definition work can separate persistent observer state from optional event-output capacity without another ownership rewrite
- keep simulator classes orchestration-focused while making their caches thinner and more obviously derived
- keep public API changes minimal unless a tiny helper falls out naturally from the internal state cleanup
- prefer a Python-owned first pass for solver state and fetched-output caches; defer richer per-work-item status modeling or a wholesale device-side solver-state ABI until it is clearly justified

## Design Constraints

- keep `InitialValueProblem` as the semantic owner of next-solve problem data only; do not let the solver-state work drift into a second IVP-like state owner
- treat persistent observer state as adjacent runtime state, but not as part of the solver-state object itself
- do not let fetched outputs (`xf`, `tf`, trajectory samples, features, event arrays) become the semantic owner of execution state simply because they are easy to cache on the host
- separate kernel build specification from runtime state; only true specialization inputs should drive program rebuilds
- treat host mirrors and device buffers as transfer or cache details of those semantic owners, not as semantic state owners in their own right
- make the internal ownership split legible in Python first, then let `_opencl/executors.py` and `_opencl/buffers.py` implement that contract
- if a per-work-item `t0` path is added in this PR, prefer the narrowest implementation that clarifies continuation semantics; do not require a full matched device-side solver-state struct unless it clearly reduces complexity

## Non-goals

- no new IVP public-surface expansion beyond incidental polish
- no observer-definition or `ObserverData` redesign in the same PR, beyond boundary-preserving internal cleanup that makes that later work easier
- no stepper-definition model work in the same PR
- no chunked trajectory streaming or storage redesign yet
- no generic SciPy interop work beyond what already landed for Python-backed IVPs
- no multi-device work
- no Random123 adoption yet
- no broad package or kernel-asset relocation

## Suggested Implementation Slices

1. Introduce explicit internal state objects or modules for solver-owned execution state and for fetched host-side output caches, and move simulator bookkeeping onto those abstractions without changing numerical behavior.
2. Refactor simulator-to-executor synchronization so IVP problem data, solver state, and cached outputs each invalidate independently and predictably.
3. Refactor executor host mirrors and common-buffer handling so they implement the new ownership contract rather than acting as parallel semantic state holders.
4. Make program rebuild and cache-key decisions follow explicit compile-time inputs rather than incidental runtime cache invalidation.
5. Verify the same ownership split across transient, feature, and trajectory paths, especially where persistent observer state and optional output buffers interact with continuation.
6. Add focused regressions for invalidation and continuation behavior before any follow-on public helper work.

## Code-Facing Checklist

- `clode/simulation/_state.py` (new internal module): add the small internal state and cache types that let `Simulator` stop owning raw `_device_*` fields directly.
- `clode/simulation/base.py`: move solver-state and fetched-output ownership onto those internal types, split IVP push/pull from solver invalidation, and make `update_x0=True` stop relying on `None` as the only state marker.
- `clode/_opencl/executors.py` and `clode/_opencl/buffers.py`: if a narrow per-work-item `t0` path proves worthwhile during implementation, thread it through the common runtime state without turning this PR into a wholesale solver-state ABI redesign.
- `clode/kernels/transient.cl`, `clode/kernels/trajectory.cl`, `clode/kernels/features.cl`, and `clode/kernels/initializeObserver.cl`: only touch these entrypoints if the PR adopts a real device-side `t0` buffer; all of them currently assume `ti = tspan[0]`.
- `clode/simulation/features.py`: distinguish build-affecting observer changes from runtime-only observer changes, and route feature invalidation through the shared state or cache layer.
- `clode/simulation/trajectory.py`: route trajectory invalidation and fetched trajectory cache ownership through the shared state or cache layer, keeping only trajectory-specific fetch and reshape logic local.
- `clode/_opencl/models.py`: extend or reuse the existing build-model types so compile-time build inputs are explicit before they become a `BuildKey`.
- `clode/_opencl/executors.py`: centralize transfer-cache invalidation, make build-trigger decisions explicit, and treat `_x0_host`, `_pars_host`, `_xf_host`, `_dt_host`, and `_tf_host` as mirrors rather than semantic state.
- `clode/_opencl/buffers.py`: keep the current kernel-facing buffer layout stable if possible, but make the logical split between problem, solver, observer, and output concerns clearer in the buffer API.
- `test/test_simulation_contracts.py`: add state-ownership and invalidation coverage.
- `test/test_opencl_models.py`, `test/test_opencl_source_builder.py`, and `test/test_opencl_runtime.py`: add compile-time-input and program-cache coverage.
- `test/core_numerics/test_transient.py`, `test/core_numerics/test_features_basicall.py`, and `test/core_numerics/test_stochastic.py`: keep continuation behavior as the numerical regression gate.

The detailed implementation plan for this PR lives in `.design/reference/solver_state_implementation_plan.md`.

## Acceptance Criteria

- split-window continuation regressions remain green for `basicall` features and seeded stochastic Euler
- there is one explicit internal home for solver/execution state distinct from IVP-owned problem data, persistent observer state, and fetched output state
- simulator-side caches are thinner and do not duplicate semantic ownership already held by the IVP or by the solver-state object
- executor host mirrors are treated as runtime transfer details rather than as parallel semantic state owners
- program rebuild triggers are easier to follow because compile-time build inputs are distinguished from runtime state changes
- invalidation rules for changing IVP data, solver parameters, `t_span`, and observer/runtime configuration are easier to follow in code and covered by focused tests
- feature and trajectory paths still behave correctly while depending on the clarified solver-state boundary
- the currently landed IVP and SciPy-callability tests remain green

## Follow-on If This Lands Cleanly

The next few high-value PRs should be separation of integration state from output and storage policy, then observer-definition cleanup around persistent state versus optional event storage, then a cleaner Python-owned stepper-definition model. Public continuation-policy helpers, richer IVP batch helpers, and broader solver interop should follow once those boundaries stop moving. A dedicated ensemble class can be revisited later only if IVP-owned batching and helper functions prove too limiting.
