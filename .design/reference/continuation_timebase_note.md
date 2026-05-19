# Continuation time-base note

Purpose: record the current live continuation semantics and the remaining state-model questions.
Read when: touching continuation helpers, attained-`tf` handoff, or solver-owned time semantics.
Update when: continuation behavior, public continuation helpers, or the internal solver-state model changes.

Historical debugging detail from the earlier observer-time rebase investigation now lives in `../archived/reference_cleanup_2026_05_13/continuation_timebase_background.md`.

## Current verified state

- The split-window continuation regressions now live in `test/core_numerics/test_features_basicall.py` and `test/core_numerics/test_stochastic.py` and pass on the stable workspace runtime.
- The live model is solver-owned absolute time. The observer no longer owns time rebasing; `finalizeObserverState` in `observer_basic_allVar.clh` is intentionally empty.
- `_opencl/buffers.py` and `_opencl/executors.py` now persist RNG state, the Box-Muller spare normal, and the prepared next-step Wiener sample across `shift_x0()` plus `set_tspan()` continuation, while clearing the prepared Wiener state on explicit `x0` or problem or solver resets.
- Exact absolute-time continuation still requires caller-managed `t_span`. For fixed-step runs, the robust handoff point is the attained `tf` returned by `get_final_time()`, not the requested endpoint.

## Remaining live design questions

- The solver now has a first-class Python-owned `SolverState` in `clode/simulation/_state.py`, but the runtime still exposes continuation through a shared requested `tspan` plus per-work-item `dt`, attained `tf`, and RNG continuation rather than a full device-side per-work-item state object.
- Per-work-item `t0` and explicit completion or error flags would make continuation and diverged-work-item bookkeeping much simpler.
- Per-work-item `t0` fits the current continuation model semantically, but the live kernels, observer initialization path, and two-pass rewind logic still assume `ti = tspan[0]`, so a real device-side `t0` is a coordinated runtime change rather than just another cache field.
- The Python-level `SolverState` first pass is now live; a full matched device-side solver-state struct remains optional future follow-through once observer-definition and stepper-definition boundaries stop moving.
- Stochastic continuation details are now preserved, but they may want a clearer per-work-item RNG-state home if Random123 or related RNG work becomes active.
- If the package wants repeated solve calls to continue absolute time more ergonomically, it likely needs an explicit continuation-policy helper rather than more implicit rebasing rules.
- Long-time fixed-step accumulation remains separate future work; fixed-step kernels still use `ti += dt`, so very long runs at small `dt` remain precision-sensitive.

## Current public-docs boundary

Public docs should stay usage-focused:

- what continues automatically
- when to use `get_final_time()`
- when to concatenate trajectory windows
- what the current limitations are
