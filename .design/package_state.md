# clODE Package State

Purpose: factual package and repo map for the live codebase.
Read when: you need to know where code belongs, which public surfaces are live, or which constraints still apply.
Update when: canonical module homes, public compatibility surfaces, packaging rules, or runtime assumptions change.

## Snapshot

- clODE is a pure-Python, PyOpenCL-only package, and runtime-critical OpenCL sources ship from `clode/kernels/` as package data.
- The public API is centered on `Simulator`, `TrajectorySimulator`, and `FeatureSimulator`, but those simulator classes should be treated as orchestration objects rather than the semantic owners of ensembles or solver state.
- Canonical public semantic homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`; `clode._opencl` is the canonical internal execution layer.
- Flat root modules such as `clode.solver`, `clode.features`, `clode.trajectory`, `clode.types`, `clode.function_converter`, `clode.xpp_parser`, and `clode.opencl_builtins` remain compatibility barrels rather than the default home for new implementation work.
- Public simulators now construct `_opencl` executors directly; the migration-era backend protocol and factory path is gone.
- `InitialValueProblem`, canonical solver-setting resolution, canonical `ObserverParams` resolution, Python-owned stepper definitions, observer definitions, and solver-state or result-cache invalidation boundaries are all live.
- Internal observer-spec resolution and `_opencl` metadata or struct helpers now consume `ObserverRuntimeSettings`, `EventOutputSettings`, `IntegrationSettings`, and `TrajectoryOutputSettings` directly; `ObserverParams` and `SolverParams` remain compatibility surfaces at public or explicitly wrapped boundaries rather than the default internal contract.
- Kernel-side observer runtime configuration now also follows that split: the shared OpenCL settings header owns `ObserverRuntimeSettings`, and internal kernel-component tests no longer need to route that contract through the public `ObserverParams` bundle.
- The runtime is explicitly single-device only, and kernels remain compile-time specialized by precision, stepper, observer, and problem shape.
- Exact shared-final-time continuation is available through `advance_tspan_to_attained_final_time()`, while diverged per-work-item final times still require explicit caller policy.
- Live numerical hardening now includes solver-owned compensated elapsed time across fixed-step and adaptive steppers, compensated `basic` and `basicall` means, inverse-linear `threshold_2` timestamps, bounded three-sample `local_max` helpers, corrected Kahan-style elapsed-time reconstruction, and fixed multi-stage stage times rebuilt from `t0 + elapsed + fractional_dt` rather than from a rounded absolute `ti`.
- The stepper wrapper boundary now preserves both failure status and accepted step width while still keeping next-step `dt` as controller state, and the runtime now surfaces per-work-item solver status, accepted-step-count, and last-accepted-step-width arrays through `get_status()`, `get_step_count()`, and `get_last_accepted_dt()`, with distinct status codes for completion, max-step exhaustion, terminal-event stop, output-capacity stop, no-progress when the requested window collapses in runtime precision or an adopted in-loop float32 step cannot advance time, and stepper failure.
- Transient, trajectory, and observer workflows are related but not interchangeable concepts: transient solves are the no-retained-output path, trajectory solves apply a retained-sample output policy, and observers remain the stateful feature or event-detection path rather than the semantic owner of all outputs.
- Public observer feature surfaces no longer report solver-owned step-count or `dt` summary diagnostics; event count remains semantically observer-owned, while some heavier observers still retain private counters needed for event geometry.
- The active planning focus is a simulation state and output ownership pass: clarify how integration settings, solver state, observer runtime settings, persistent observer state, trajectory output policy, and fetched outputs are modeled on the Python side versus the `_opencl` execution side before taking on more observer or output-surface growth.

## Contributor Routing

| Area | Canonical homes | Notes |
| --- | --- | --- |
| Problem authoring and RHS ingestion | `clode/problem/*` | The public problem API centers `InitialValueProblem`; lower-level support types stay internal under `_core`. |
| Simulation and orchestration | `clode/simulation/*` | New implementation work belongs here even when compatibility barrels re-export it. |
| Observers and feature metadata | `clode/observers/*`, `clode/kernels/observers/*.clh` | Python observer definitions drive feature names and build metadata; runtime-specific layouts live in `_opencl/observer_metadata.py`. |
| Runtime, build, and dispatch | `clode/_opencl/*`, `clode/kernels/*` | `_opencl` owns runtime, build, buffer, cache, and dispatch mechanics; treat `clode/kernels/odedriver.cl` as deferred design context, not active runtime code. |
| Public runtime selection and logging | `clode/runtime/*` | Runtime selection is explicit single-device selection plus stdlib logging helpers. |
| Tests and evidence | `test/core_numerics/`, `test/kernel_components/`, `test/test_simulation_contracts.py`, `tools/run_test_bundle.py` | Prefer component tests for helper or build contracts and numerics for exact-solution evidence. |
| Public docs and repo surface | `docs/`, `README.md`, `paper/paper.md` | Keep user-facing wording current and non-historical. |

## Live Constraints And Debt

- The runtime is single-device only. Any future multi-device execution needs a dedicated API and execution model.
- Simulators still carry compatibility delegates and cached mirrors alongside the IVP; richer IVP-side batch helpers such as grids, random sampling, and quasi-random sampling are still missing.
- `SolverState` is a useful first pass, but there is still no matched device-side per-work-item solver-state object for current time or richer stepping diagnostics beyond the surfaced per-item status, step-count, and last-accepted-step-width buffers.
- `advance_tspan_to_attained_final_time()` only covers the representable shared-final-time case. `shift_tspan()` remains the requested-window continuation tool, not the exact attained-time path.
- Public `SolverParams` still mixes integration policy with trajectory-output policy, `ObserverParams` still mixes runtime thresholds with event-output capacity at the compatibility surface, and `FeatureSimulator` still exposes legacy `observer_*` inputs alongside canonical `ObserverParams`.
- The flat root modules plus `SolverParams` and `ObserverParams` are now explicitly treated in code as compatibility surfaces rather than as primary semantic owners, but the eventual public API audit and cleanup is still deferred until the narrower owner model settles further.
- There is still no explicit Python-side owner split between persistent observer state, observer runtime settings, event-output policy, and fetched feature or trajectory outputs; that boundary is currently spread across simulator caches, `_opencl` metadata, and transfer-cache helpers.
- The public solver-diagnostics API is still intentionally fine-grained. A bundled stats object is deferred until the remaining fields and work-metric semantics settle.
- Some observer structs still carry private step-count bookkeeping for event geometry or running means, but public observer outputs no longer expose solver-owned step-count or `dt` summary diagnostics; event counts remain semantically observer-owned.
- Heavier observer-state footprint and register-pressure work remains backlog rather than active scope.
- Numerical validation and public evidence still need broader exact-solution and convergence coverage beyond the landed stable-linear transient evidence slices for RK4 global-error convergence and Dormand-Prince tolerance refinement to match the current claims.
- RNG continuation details still live in separate buffers rather than a clearer named per-work-item state model.
- Packaging, citation, and other repo-surface cleanup remain lower priority than numerical, runtime, and UX work.

## Deep Dives

- Layout and semantic ownership: `.design/reference/semantic_layout_audit.md`
- Solver-state boundary and deferred follow-through: `.design/reference/solver_state_implementation_plan.md`
- Continuation and time-base behavior: `.design/reference/continuation_timebase_note.md`
- Float32 numerics and helper guardrails: `.design/reference/single_precision_numerics_note.md`
- Testing strategy and bundle intent: `.design/reference/testing_audit.md`
- PyOpenCL leverage and runtime helper opportunities: `.design/reference/pyopencl_leverage_audit.md`
- Public docs and publication surfaces: `.design/reference/docs_layout_plan.md`, `.design/reference/public_surfaces_plan.md`, `.design/reference/joss_audit.md`
