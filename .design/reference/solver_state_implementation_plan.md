# Solver-State Implementation Plan

Purpose: detailed code-facing implementation map for the active solver-state ownership PR.
Read when: implementing `.design/next_pr.md`, breaking the work into commits, or deciding which modules and tests belong in scope.
Update when: the file targets, phase ordering, or accepted internal shapes for the solver-state PR change.

## Bottom line

- Keep the public API stable while replacing the current overlapping ownership model with explicit internal state and cache boundaries.
- Reuse the existing `_opencl` build-key and program-cache pipeline instead of inventing a second build system; make compile-time inputs explicit before they become a `BuildKey`.
- Prefer semantic cleanup first and kernel-ABI churn second. If the current buffer layout can survive this PR, keep it.
- Treat this PR as the contract-setting pass that the next three workstreams need: output-policy separation, observer-definition cleanup, and stepper-definition cleanup.

## Current ownership hot spots

### `clode/simulation/base.py`

- Owns IVP orchestration and also owns `_device_initial_state`, `_device_parameters`, `_device_final_state`, `_device_dt`, `_device_tf`, `_t_span`, and `_cl_program_is_valid`.
- `_sync_problem_data_from_ivp()` and `_sync_ivp_from_device_problem_data()` currently mix IVP ownership with runtime cache refresh and device pullback.
- `transient(update_x0=True)` advances device-side `x0`, then uses `None`-based cache invalidation rather than one explicit state boundary.

### `clode/simulation/features.py`

- Owns feature fetch caching plus observer invalidation and rebuild decisions.
- `set_observer()` and `set_observer_parameters()` currently combine build invalidation, observer-state invalidation, and feature-cache invalidation.

### `clode/simulation/trajectory.py`

- Owns trajectory fetch caches and reshaping, but still relies on base-class transient cache ownership and device-state invalidation.

### `clode/_opencl/executors.py`

- Owns `_tspan`, `_solver_params`, `_program_bundle`, `_buffers`, `_trajectory_buffers`, `_feature_buffers`, and a parallel host-mirror layer (`_x0_host`, `_pars_host`, `_xf_host`, `_dt_host`, `_tf_host`, `_rng_state_host`).
- `set_problem_data()`, `set_solver_params()`, `set_tspan()`, `set_observer()`, and `set_observer_params()` each invalidate different overlapping subsets of those mirrors and buffers.
- Feature and trajectory subclasses add yet another layer of fetch-cache invalidation and buffer lifecycle rules.

### `clode/_opencl/buffers.py`

- `CommonBuffers` still mixes problem data, solver state, RNG state, and transient outputs into one buffer bundle.
- That layout is fine as a kernel ABI, but it is not yet reflected as clear semantic groups in Python.

### `_opencl` build pipeline

- `_opencl/models.py`, `_opencl/source_builder.py`, and `_opencl/program_cache.py` already hold the real compile-time artifacts.
- The missing piece is an explicit, code-facing distinction between build-affecting inputs and runtime-only state changes at the simulator and executor boundary.

## Recommended internal shapes

These names are suggestions, not frozen API.

### `clode/simulation/_state.py` (new internal module)

- `SolverState`: semantic execution state for the current requested solve window and its continuation-relevant per-work-item fields.
- `TransientCache`: fetched host-side final-state, `dt`, and `tf` caches only.
- `TrajectoryCache`: fetched host-side trajectory arrays and `n_stored` only.
- `FeatureCache`: fetched feature table and derived feature-count cache only.

Recommended rule:

- IVP still owns problem definition and next-solve input data.
- Solver state owns requested `t_span`, semantic current `dt`, attained `tf`, RNG continuation, and any future status flags.
- Output caches own fetched host arrays only.

### `clode/_opencl/models.py`

- Add a thin `BuildInputs` or `KernelBuildInputs` value object if the simulator and executor layers need something explicit before a `BuildKey` exists.
- Keep `BuildKey` as the backend-ready cache key used by `SourceBuilder` and `ProgramCache`.

Recommended rule:

- If a change does not affect the eventual `BuildKey`, it should not force a program rebuild.

## Phase plan

### Phase 1: Introduce internal state and cache shells

Primary files:

- `clode/simulation/_state.py`
- `clode/simulation/base.py`
- `clode/simulation/features.py`
- `clode/simulation/trajectory.py`

Tasks:

- Add the new internal state or cache dataclasses.
- Move simulator-side `_device_final_state`, `_device_dt`, `_device_tf`, feature caches, and trajectory caches onto those objects.
- Keep behavior unchanged; this phase is about ownership and naming, not semantics.

Exit criteria:

- `Simulator` no longer uses raw `_device_*` fields as the primary semantic boundary for fetched outputs.
- Feature and trajectory fetch caches are visibly separate from solver-state ownership.

### Phase 2: Make the IVP push/pull boundary explicit

Primary files:

- `clode/simulation/base.py`
- optionally `clode/simulation/_state.py`

Tasks:

- Replace the current implicit sync flow with explicit helpers for pushing IVP data to the executor and pulling device-updated problem data back into the IVP when needed.
- Make `update_x0=True` stop relying on `self._device_initial_state = None` as the only signal that the IVP and runtime have diverged.
- Decide one narrow policy for this PR:
  either sync the IVP immediately after `shift_x0()` or track one explicit “problem data needs pullback” flag.

Recommended choice:

- Prefer one explicit pending-pull flag over unconditional eager download if eager sync would force unnecessary readback.

Exit criteria:

- There is one obvious place in the code where IVP ownership yields to runtime-updated initial state and one obvious place where it is pulled back.

### Phase 3: Separate build invalidation from runtime invalidation

Primary files:

- `clode/_opencl/models.py`
- `clode/_opencl/source_builder.py`
- `clode/_opencl/executors.py`
- `clode/simulation/base.py`
- `clode/simulation/features.py`

Tasks:

- Surface the compile-time inputs explicitly before build time.
- Make simulator or executor rebuild decisions follow those inputs rather than `_cl_program_is_valid` toggles alone.
- Keep runtime-only changes out of that path: `t_span`, `dt`, threshold values, and fetched outputs should not look like build-state changes.

Key cases to preserve:

- `set_observer()` should rebuild features kernels.
- `set_observer_parameters(max_event_timestamps=...)` should rebuild features kernels.
- threshold or neighborhood-radius changes should reset observer runtime state but not trigger a rebuild.
- `set_tspan()` and `set_solver_parameters()` should not rebuild transient kernels.

Exit criteria:

- The code can answer “why did this rebuild?” in one place.
- The code can answer “why did this only clear runtime state?” in one place.

### Phase 4: Clean executor transfer-cache ownership

Primary files:

- `clode/_opencl/executors.py`
- `clode/_opencl/buffers.py`

Tasks:

- Group executor invalidation helpers by concern instead of by ad hoc field resets.
- Treat `_x0_host`, `_pars_host`, `_xf_host`, `_dt_host`, `_tf_host`, and related arrays as mirrors only.
- Preserve the current `CommonBuffers` layout if possible, but make the logical groups explicit in Python:
  - problem data
  - solver state
  - RNG continuation state
  - transient outputs
  - trajectory outputs
  - feature outputs and observer state

Recommended choice:

- Avoid a kernel-argument reorder or buffer-ABI redesign in this PR unless it falls out almost for free.

Exit criteria:

- Executor setters clear only the mirrors and buffers they actually invalidate.
- Buffer allocation and mirror invalidation are easier to follow than the current overlapping resets.

### Phase 5: Reconcile feature and trajectory subclasses with the new state boundary

Primary files:

- `clode/simulation/features.py`
- `clode/simulation/trajectory.py`
- `clode/_opencl/executors.py`

Tasks:

- Make feature execution depend on explicit observer-state invalidation rather than on broad simulator cache clearing.
- Make trajectory execution depend on output-policy and trajectory-buffer invalidation rather than on broad transient cache semantics.
- Keep observer initialization state tied to observer-state and runtime mutations, not to unrelated transient cache changes.

Exit criteria:

- Feature and trajectory code paths still behave like thin orchestration layers over shared base semantics.
- Subclass-local caches are only the caches that genuinely belong to that subclass.

### Phase 6: Focused regression coverage

Primary test files:

- `test/test_simulation_contracts.py`
- `test/test_opencl_models.py`
- `test/test_opencl_source_builder.py`
- `test/test_opencl_runtime.py`
- `test/core_numerics/test_transient.py`
- `test/core_numerics/test_features_basicall.py`
- `test/core_numerics/test_stochastic.py`

Add or expand tests for:

- changing solver parameters resets runtime state without causing rebuilds
- changing `t_span` invalidates solver-state outputs without invalidating IVP-owned problem data
- `update_x0=True` yields a coherent next-solve initial-state story
- observer changes split correctly into rebuild-affecting versus runtime-only mutations
- program-cache reuse still keys off build inputs rather than runtime-only state changes
- split-window continuation remains green for transient and feature paths

## Commit slicing recommendation

1. Add internal state or cache dataclasses and move simulator fetch caches onto them.
2. Make IVP push or pull behavior explicit.
3. Introduce explicit build-input handling and rebuild rules.
4. Clean executor mirror invalidation and buffer grouping.
5. Reconcile feature and trajectory paths with the new state model.
6. Add focused tests and then rerun the continuation regression slice.

## Explicit non-goals for this plan

- Do not split `SolverParams` into integration and output-policy objects in this PR; that is the next PR.
- Do not redesign observer definitions in this PR; only leave behind a cleaner boundary for that work.
- Do not broaden this PR into generic SciPy interop, new batch-helper APIs, or multi-device abstractions.
- Do not add a backend-agnostic execution layer.

## Open implementation choices to resolve early

### Status flags

- Keep the first pass narrow. If true per-work-item status flags require broader kernel changes than expected, prefer a minimal internal status model over a rich new error taxonomy.

### IVP sync after `shift_x0()`

- Prefer an explicit pending-pull flag if eager synchronization would add unnecessary host readback.

### Build-spec type

- Prefer extending `_opencl/models.py` with a thin precursor to `BuildKey` rather than inventing a second build-spec structure elsewhere.
