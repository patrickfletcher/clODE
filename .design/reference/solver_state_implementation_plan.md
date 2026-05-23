# Solver-State Implementation Plan

Purpose: record what landed in the first-pass solver-state ownership cleanup and what it intentionally deferred.
Read when: follow-on PRs need the exact current boundary between IVP-owned problem data, solver state, persistent observer state, and fetched outputs.
Update when: the solver-state boundary changes materially or this note becomes historical enough to archive.

## Bottom line

- The first-pass solver-state cleanup landed.
- IVP owns next-solve problem data.
- `clode/simulation/_state.py` now owns Python-side solver state plus fetched-output caches.
- `_opencl/executors.py` now treats host mirrors as transfer caches rather than semantic state owners.
- The stepper wrapper boundary now preserves failure status and accepted step width explicitly, but no matched device-side per-work-item solver-state struct was added in this pass.

## Landed internal shape

### `clode/simulation/_state.py`

- `SolverState` owns the requested window, current continued time, current `dt`, attained `tf`, and the minimal sync marker needed for runtime-updated `x0` pullback.
- `TransientCache`, `FeatureCache`, and `TrajectoryCache` own fetched host arrays only.
- Feature and trajectory caches now carry an explicit `has_result` marker so stale device buffers are not silently re-downloaded after invalidation.

### `clode/simulation/base.py`

- IVP remains the sole semantic owner of next-solve problem data on the Python side.
- Simulator mutators now invalidate solver-owned results and subclass caches explicitly instead of relying on raw `_device_*` mirrors.
- `update_x0=True` now continues solver-owned time and marks IVP problem data for on-demand pullback rather than using `None`-only cache conventions.

### `clode/_opencl/executors.py` and `clode/_opencl/buffers.py`

- Executor-side host mirrors are grouped as transfer caches by concern rather than as ad hoc fields.
- Changing problem data or solver parameters resets per-item `dt` to the requested solver `dt`.
- `tspan`-only changes preserve continuation `dt` while invalidating prior `xf` and `tf` fetches.
- Trajectory and feature subclasses now keep their own transfer caches without retaining stale retired-buffer lists.

## Verified coverage

- `test/test_simulation_contracts.py` covers on-demand IVP sync plus stale-result invalidation.
- `test/test_opencl_executors.py` covers executor-side `dt` reset versus continuation-preserve semantics.
- `test/core_numerics/test_transient.py`, `test/core_numerics/test_features_basicall.py`, and `test/core_numerics/test_stochastic.py` remain the continuation regression gate.

## Intentionally deferred

- device-side per-work-item `t0`
- richer per-work-item completion/error status
- a matched device-side solver-state struct
- solver-owned device buffers for step counts, accepted step width, and other per-item stepping diagnostics
- separation of integration state from output/storage policy (`max_store`, `nout`, event capacity)
- continuation-state guardrails for diverged per-work-item final times and any later public continuation helper
- migration of legacy step/time diagnostic ownership out of observer structs and outputs
- observer-definition cleanup beyond preserving the boundary for later work, which has now landed separately

## Handoff to the next PR

- Treat the landed solver-state boundary as stable enough to build on, not as something to reopen wholesale.
- The next structural cleanup is to make solver-owned per-work-item stepping state explicit on the device and runtime sides: step counts, accepted or current `dt`, time-base values, and failure status.
- Strip matching legacy step or time diagnostics out of observer ownership while keeping observer-specific interpolation and event geometry intact.
- Keep the current owner split intact in follow-on work:
  - IVP owns next-solve problem data
  - solver state owns execution progress, accepted-step and status facts, and continuation semantics
  - persistent observer state stays adjacent but separate
  - fetched outputs and transfer caches remain derived
