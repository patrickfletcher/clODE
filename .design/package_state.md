# clODE Package State

Purpose: factual package and repo map for the live codebase.
Read when: you need to know where code belongs, which public surfaces are live, or which constraints still apply.
Update when: canonical module homes, public compatibility surfaces, packaging rules, or runtime assumptions change.

## Snapshot

- clODE is a pure-Python, PyOpenCL-only package, and runtime-critical OpenCL sources ship from `clode/kernels/` as package data.
- The public API is centered on `Simulator`, `TrajectorySimulator`, and `FeatureSimulator`, but those simulator classes should be treated as orchestration objects rather than the semantic owners of ensembles or solver state.
- Canonical public semantic homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`; `clode._opencl` is the canonical internal execution layer.
- Flat root modules such as `clode.solver`, `clode.features`, `clode.trajectory`, `clode.types`, `clode.function_converter`, `clode.xpp_parser`, and `clode.opencl_builtins` remain compatibility barrels rather than the default home for new implementation work.
- `clode/function_converter.py` and `clode/xpp_parser.py` are thin compatibility re-export shims; the maintained implementations now live in `clode/problem/python.py` and `clode/problem/xpp.py`.
- Public simulators now construct `_opencl` executors directly; the migration-era backend protocol and factory path is gone.
- `InitialValueProblem`, canonical solver-setting resolution, canonical `ObserverParams` resolution, Python-owned stepper definitions, observer definitions, and solver-state or result-cache invalidation boundaries are all live.
- Internal observer-spec resolution and `_opencl` metadata or struct helpers now consume `ObserverRuntimeSettings`, `EventOutputSettings`, `IntegrationSettings`, and `TrajectoryOutputSettings` directly; `ObserverParams` and `SolverParams` remain compatibility surfaces at public or explicitly wrapped boundaries rather than the default internal contract.
- Kernel-side observer runtime configuration now also follows that split: the shared OpenCL settings header owns `ObserverRuntimeSettings`, and internal kernel-component tests no longer need to route that contract through the public `ObserverParams` bundle.
- Current naming direction is explicit: `SolverState` is solver-owned live execution state, `TrajectoryOutput` is fetched retained-sample data, device-side `ObserverState` is persistent observer-private runtime state, and public `ObserverOutput` is the fetched feature or event readout produced from that state rather than the persistent state object itself.
- The runtime is explicitly single-device only, and kernels remain compile-time specialized by precision, stepper, observer, and problem shape.
- Exact shared-final-time continuation is available through `advance_tspan_to_attained_final_time()`, while diverged per-work-item final times still require explicit caller policy.
- Live numerical hardening now includes solver-owned compensated elapsed time across fixed-step and adaptive steppers, compensated summary means, compensated trajectory and auxiliary means across the semantic event observers, compensated Schmitt downstate means for `active dip`, inverse-linear `normalized_schmitt_trigger` timestamps, bounded three-sample `local_max` helpers, corrected Kahan-style elapsed-time reconstruction, and fixed multi-stage stage times rebuilt from `t0 + elapsed + fractional_dt` rather than from a rounded absolute `ti`.
- The stepper wrapper boundary now preserves both failure status and accepted step width while still keeping next-step `dt` as controller state, and the runtime now surfaces per-work-item solver status, accepted-step-count, and last-accepted-step-width arrays through `get_status()`, `get_step_count()`, and `get_last_accepted_dt()`, with distinct status codes for completion, max-step exhaustion, terminal-event stop, output-capacity stop, no-progress when the requested window collapses in runtime precision or an adopted in-loop float32 step cannot advance time, and stepper failure.
- Transient, trajectory, and observer workflows are related but not interchangeable concepts: transient solves are the no-retained-output path, trajectory solves apply a retained-sample output policy, and observers remain the stateful feature or event-detection path rather than the semantic owner of all outputs.
- Public observer feature surfaces no longer report solver-owned step-count or `dt` summary diagnostics; event count remains semantically observer-owned, while some heavier observers still retain private counters needed for event geometry.
- The summary-observer family proof slice is now landed: `Observer.summary` plus `SummaryObserverSelection` specialize summary feature schemas and persistent state layouts to the requested subset. Event observers include `Observer.threshold_crossing`, `Observer.normalized_threshold_crossing`, `Observer.schmitt_trigger`, `Observer.normalized_schmitt_trigger`, `Observer.local_max`, and `Observer.normalized_neighborhood_return`. `ThresholdCrossingConfig`, `SchmittTriggerConfig`, `LocalMaximumConfig`, and `NeighborhoodReturnConfig` are live semantic config surfaces; all current event-triggering families now expose the canonical `min_amp` and `max_event_count` controls, threshold, Schmitt, and neighborhood configs expose both `event_var` and `feature_var` where their kernels split trigger and measurement channels, and `LocalMaximumConfig` remains intentionally one-channel. `ObserverParams` and legacy `observer_*` keywords remain compatibility adapters. Shared accepted-step history update helpers (`advanceAcceptedStepHistory2`, `advanceAcceptedStepHistory3`, and their `ByVariable` variants) are now live in `observers.cl` and consumed by the observer families.

## Contributor Routing

| Area | Canonical homes | Notes |
| --- | --- | --- |
| Problem authoring and RHS ingestion | `clode/problem/*` | The public problem API centers `InitialValueProblem`; lower-level support types stay internal under `_core`. |
| Simulation and orchestration | `clode/simulation/*` | New implementation work belongs here even when compatibility barrels re-export it. |
| Observers and feature metadata | `clode/observers/*`, `clode/kernels/observers/*.clh` | Python observer definitions drive feature names and build metadata; runtime-specific layouts live in `_opencl/observer_metadata.py`. |
| Runtime, build, and dispatch | `clode/_opencl/*`, `clode/kernels/*` | `_opencl` owns runtime, build, buffer, cache, and dispatch mechanics. |
| Public runtime selection and logging | `clode/runtime/*` | Runtime selection is explicit single-device selection plus stdlib logging helpers. |
| Tests and evidence | `test/core_numerics/`, `test/kernel_components/`, `test/test_simulation_contracts.py`, `tools/run_test_bundle.py` | Prefer component tests for helper or build contracts and numerics for exact-solution evidence. |
| Public docs and repo surface | `docs/`, `README.md`, `paper/paper.md` | Keep user-facing wording current and non-historical. |

## Live Constraints And Debt

- The runtime is single-device only. Any future multi-device execution needs a dedicated API and execution model.
- Simulators still carry compatibility delegates and cached mirrors alongside the IVP; richer IVP-side batch helpers such as grids, random sampling, and quasi-random sampling are still missing.
- `SolverState` is a useful first pass, but there is still no matched device-side per-work-item solver-state object for current time or richer stepping diagnostics beyond the surfaced per-item status, step-count, and last-accepted-step-width buffers.
- `advance_tspan_to_attained_final_time()` only covers the representable shared-final-time case. `shift_tspan()` remains the requested-window continuation tool, not the exact attained-time path.
- Public `SolverParams` still mixes integration policy with trajectory-output policy, `ObserverParams` still mixes canonical event-trigger controls, runtime thresholds, and event-output capacity at the compatibility surface, and `FeatureSimulator` still exposes legacy `observer_*` inputs alongside canonical `ObserverParams`. Current observer UX now makes that debt concrete: `max_event_count` and `min_amp` are family-wide semantic controls across the event-triggering observers, but compatibility-only fields such as `eps_dx` and `min_imi` still sit beside those canonical controls inside the same broad bundle.
- The flat root modules plus `SolverParams` and `ObserverParams` are now explicitly treated in code as compatibility surfaces rather than as primary semantic owners, but the eventual public API audit and cleanup is still deferred until the narrower owner model settles further.
- The owner split is clearer than it was before the last refactor, and the summary family now has a first reusable declaration path, but there is still no observer-wide Python-side model that packages event condition, persistent observer state, runtime settings, event-output policy, and fetched observer outputs into one easy-to-extend definition layer across the heavier event families.
- Event-observer selection granularity is still coarse, but the one-pass summary family now supports explicit state, auxiliary, and slope reduction subsets through `SummaryObserverSelection`; those summary subsets are currently build-specialized rather than hot-swapped through narrower runtime settings.
- Whether to expose raw observer state directly, and whether `ObserverOutput` should remain the long-lived public name for the current feature or event readout object, is intentionally deferred until the observer authoring and packaging model is settled enough to compare alternatives without compatibility pressure.
- The public solver-diagnostics API is still intentionally fine-grained. A bundled stats object is deferred until the remaining fields and work-metric semantics settle.
- Some observer structs still carry private step-count bookkeeping for event geometry or running means, but public observer outputs no longer expose solver-owned step-count or `dt` summary diagnostics; event counts remain semantically observer-owned.
- **Observer readout architecture**: canonical event observers expose event geometry plus observer-local summary/readout groups in observer-specific schemas. `local_max` remains the public extrema family with dual max/min event streams and IMI/amplitude outputs. `max_event_count` is already the general event-loop limiter across the current event-triggering families. All current event-triggering semantic configs now expose `min_amp` as an event-var range gate. Canonical threshold and Schmitt configs both split trigger geometry on `event_var` from extrema/amplitude measurement on `feature_var`; threshold keeps its single-stream trigger semantics, while Schmitt adds duration/duty and family-local `active dip` readouts. `normalized_neighborhood_return` likewise exposes `feature_var` on its semantic config and uses it for maxima/amplitude tracking while `event_var` still owns anchor geometry and the event-var gate. Kernel-side helper hardening is now standardized across the semantic event observers (shared K=3 history updates, compensated trajectory/auxiliary means, and family-consistent interpolation primitives). The remaining open question is readout bundle/selectability and schema packaging, not baseline time/mean helper alignment.
- Heavier observer-state footprint and register-pressure work remains backlog rather than active scope.
- Numerical validation and public evidence still need broader exact-solution and convergence coverage beyond the landed stable-linear transient evidence slices for RK4 global-error convergence and Dormand-Prince tolerance refinement to match the current claims.
- The XPP ingestion path remains a lightweight line-oriented parser/rewriter (`clode/problem/xpp.py`) rather than a grammar-backed parser, so edge-case syntax handling is still an explicit hardening area.
- Kernel-component evidence exists and is useful, but it is still concentrated in `test/kernel_components/test_kernel_math.py`; broader component-level coverage remains backlog work.
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
- Observer event and feature variable semantics: `.design/reference/observer_event_feature_variable_contract.md`
