# Continuation time-base note

Purpose: record the current live continuation semantics and the remaining state-model questions.
Read when: touching continuation helpers, attained-`tf` handoff, or solver-owned time semantics.
Update when: continuation behavior, divergent-time guardrails, or the internal solver-state model changes.

Historical debugging detail from the earlier observer-time rebase investigation now lives in `../archived/reference_cleanup_2026_05_13/continuation_timebase_background.md`.

## Current verified state

- The split-window continuation regressions now live in `test/core_numerics/test_features_summary.py` and `test/core_numerics/test_stochastic.py` and pass on the stable workspace runtime.
- The live model is solver-owned absolute time. The observer no longer owns time rebasing; `finalizeObserverState` in `observer_summary.clh` is intentionally empty.
- `_opencl/buffers.py` and `_opencl/executors.py` now persist RNG state, the Box-Muller spare normal, and the prepared next-step Wiener sample across `shift_x0()` plus `set_tspan()` continuation, while clearing the prepared Wiener state on explicit `x0` or problem or solver resets.
- Exact absolute-time continuation still requires caller-managed `t_span`. For fixed-step runs, the robust handoff point is the attained `tf` returned by `get_final_time()`, not the requested endpoint.
- Exact shared-window continuation is only representable when the ensemble agrees on one attained final time. Once work-items diverge in `tf`, there is no single correct next shared `t_span` update.
- The live kernels now keep solve-local compensated elapsed time and observer-relative elapsed bookkeeping internally for both fixed-step and adaptive solves, fixed multi-stage methods now rebuild their internal stage times from the same `t0 + elapsed + fractional_dt` helper path as adaptive steppers, and the stepper wrapper boundary now preserves accepted step width plus failure status explicitly.
- That stepping state is still reconstructed from the requested `tspan[0]` on each solve rather than being exposed as a new persisted host-visible current-time or per-item status buffer.

## Remaining live design questions

- The solver now has a first-class Python-owned `SolverState` in `clode/simulation/_state.py`, but the runtime still exposes continuation through a shared requested `tspan` plus per-work-item `dt`, attained `tf`, and RNG continuation rather than a full device-side per-work-item state object.
- The simulator layer now has a narrow helper for the representable case: `advance_tspan_to_attained_final_time()` proves one shared attained final time before updating the shared requested window. A broader built-in divergent-time continuation policy should still stay deferred because, with only a shared requested `tspan`, any such policy would be an approximation or arbitrary choice.
- Per-work-item `t0` and explicit completion or error flags would make continuation and diverged-work-item bookkeeping much simpler.
- Per-work-item `t0` fits the current continuation model semantically, but the live kernels, observer initialization path, and two-pass rewind logic still assume `ti = tspan[0]`, so a real device-side `t0` is a coordinated runtime change rather than just another cache field.
- The Python-level `SolverState` first pass is now live; a full matched device-side solver-state struct remains the next logical follow-through so per-item step counts, accepted step widths, and failure status stop being implicit or observer-adjacent.
- Stochastic continuation details are now preserved, but they may want a clearer per-work-item RNG-state home if Random123 or related RNG work becomes active.
- If the package wants repeated solve calls to continue exact absolute time after diverged per-item final times, it likely needs a per-work-item current-time or `t0` model rather than another shared-window helper.
- Long-time time-base work remains separate future work; the live fixed-step and adaptive paths both reconstruct time from `tspan[0]` plus a compensated solve-relative elapsed pair. Very long runs at small `dt` are still limited by float32 absolute-time spacing when an output must be stored as one absolute timestamp.
- The current narrow implementation slice now keeps one live time-base story explicit: compensated solve-relative elapsed time across steppers, while a persisted per-work-item current-time or `t0` model remains deferred until exact diverged-time continuation becomes a stronger requirement.
- None of the current time/interpolation helper candidates requires the general nonlinear system solve machinery that later implicit methods will need; those larger dependencies should remain deferred with the implicit-stepper work itself.

## Current public-docs boundary

Public docs should stay usage-focused:

- what continues automatically
- when `advance_tspan_to_attained_final_time()` is the right tool
- when to use `get_final_time()`
- when to concatenate trajectory windows
- what the current limitations are
- what compensated time bookkeeping and structured comparison baselines do and do not fix in float32
