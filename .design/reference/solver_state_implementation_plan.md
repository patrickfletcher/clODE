# Solver-State Implementation Plan

Purpose: record what landed in the first-pass solver-state ownership cleanup and what it intentionally deferred.
Read when: follow-on PRs need the exact current boundary between IVP-owned problem data, solver state, persistent observer state, and fetched outputs.
Update when: the solver-state boundary changes materially or this note becomes historical enough to archive.

## Bottom line

- The first-pass solver-state cleanup landed.
- Follow-on slices have since landed solver-owned per-work-item status, accepted-step-count, and last-accepted-step-width buffers plus public `get_status()`, `get_step_count()`, and `get_last_accepted_dt()` fetch paths, and now surface a maintained `NO_PROGRESS` status both when a positive requested window collapses to zero in runtime precision and for the currently adopted in-loop float32 time-stall cases.
- IVP owns next-solve problem data.
- `clode/simulation/_state.py` now owns Python-side solver state plus fetched-output caches.
- `_opencl/executors.py` now treats host mirrors as transfer caches rather than semantic state owners.
- The stepper wrapper boundary now preserves failure status and accepted step width explicitly, but there is still no matched device-side per-work-item solver-state struct for the remaining time-base diagnostics or later work metrics.

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
- a matched device-side solver-state struct for current-time or other remaining time-base values
- solver-owned device buffers or a matched state object for work metrics beyond status, step counts, and last accepted step width
- separation of integration state from output/storage policy (`max_store`, `nout`, event capacity)
- continuation-state guardrails for diverged per-work-item final times and any later public continuation helper
- a read-only bundled solver stats surface and later work metrics such as RHS evaluation counts, Jacobian evaluations, or linear-solver work once those quantities are semantically stable across steppers
- deeper observer-state cleanup beyond retiring public observer-owned step-count and `dt` summary outputs

## Stable follow-on guardrails

- Use the root `.design` docs for the active PR sequence. This note is the solver-state boundary record, not the current task tracker.
- Treat the landed solver-state boundary as stable enough to build on, not as something to reopen wholesale.
- The deterministic solver-owned `NO_PROGRESS` policy for the currently adopted in-loop time-stall cases is now covered by maintained public-contract tests, so any later follow-through should revisit deeper current-time modeling only if a concrete continuation or work-metric need justifies it.
- Keep observer public outputs free of solver-owned step or time diagnostics while leaving observer-specific interpolation geometry and any private sample counters intact.
- Keep the public API fine-grained for now. A bundled stats object should wait until the remaining fields and any SciPy-style work metrics have one coherent cross-stepper definition.
- Keep the current owner split intact in follow-on work:
  - IVP owns next-solve problem data
  - solver state owns execution progress, step counts, accepted-step and status facts, and continuation semantics
  - persistent observer state stays adjacent but separate
  - fetched outputs and transfer caches remain derived
