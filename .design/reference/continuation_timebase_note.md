# Continuation time-base note

Purpose: record the live continuation semantics and the solver or observer time-ownership contract.
Read when: touching continuation helpers, per-item solver time state, or observer elapsed bookkeeping.
Update when: continuation behavior, public simulator semantics, or observer time ownership changes.

Historical debugging detail from the earlier observer-time rebase investigation now lives in `../archived/reference_cleanup_2026_05_13/continuation_timebase_background.md`.

## Bottom line

- Absolute simulation time is the primary semantic time concept. Chunk-local elapsed time remains solver-private numerics.
- The runtime now persists a per-work-item solver-owned `t0` buffer. Kernels start each solve from `t0[i]` and use the shared `tspan[1] - tspan[0]` only as the requested chunk duration.
- `set_tspan()` is the explicit reset path: it updates the nominal requested window and resets device `t0` to the requested start.
- `shift_tspan()` is the continuation path: it advances the nominal requested window by one nominal duration and copies device `tf -> t0` for the next solve.
- Public `transient()`, `trajectory()`, and `features()` calls with `update_x0=True` now continue both state and time by default by performing the `xf -> x0` and `tf -> t0` handoff together.
- `advance_tspan_to_attained_final_time()` is removed. Exact attained-time continuation now lives in the `shift_tspan()` plus per-item `t0` model instead of a separate shared-final-time helper.
- Summary and semantic event observers no longer re-sum accepted-step elapsed time on every accepted step. They consume solver-owned chunk elapsed directly and persist only chunk-start elapsed offsets plus event-local time state.
- Monotone time bookkeeping now uses the Kahan-style `compensatedTime*` helpers. Means and other arbitrary sums stay on the Neumaier or integral helper path.

## Fast path

- Read `## Bottom line`, `## Live model`, `## Settled public semantics`, and `## Evidence`.
- Use `.design/reference/single_precision_numerics_note.md` for the Kahan-versus-Neumaier guardrail.
- Use `.design/reference/chunked_execution_audit.md` only for follow-on chunk-loop work now that the time-base contract is landed.

## Live model

### Solver-owned time state

| Concept | Owner | Persistence | Meaning |
| --- | --- | --- | --- |
| `t0[i]` | solver state | per work item, across continued chunks | absolute time at the start of the current solve chunk |
| `solveElapsed`, `solveElapsedCorrection` | solver local | one solve call | chunk-local compensated elapsed time |
| `ti` | solver local | current accepted step | absolute accepted time reconstructed from `t0[i]` plus the chunk-local compensated elapsed pair |
| `tf[i]` | solver output/state | per work item after each solve | attained absolute final time for the just-completed solve |
| `tspan` | solver input | per solve call | nominal requested window; kernels use its duration, host APIs expose it as the requested window |

### Observer-owned time state

| Concept | Owner | Persistence | Meaning |
| --- | --- | --- | --- |
| `tStart` | observer state | across chunks until observer reinitialization | absolute start time of the observation window |
| `elapsedAtChunkStart`, `elapsedAtChunkStartCorrection` | observer state | across chunks | observation-window elapsed time at the start of the current chunk |
| `elapsedbuffer[3]` | observer state | rolling K=3 history | observation-window elapsed coordinates for the last three accepted samples |
| `elapsedLastEvent` | observer state | across chunks | observation-window elapsed time of the last event |
| `elapsedThisDown` | Schmitt families only | across chunks | observation-window elapsed time of the last down transition |

### Consequences

- `tStart` remains the absolute start of the observation window, not the current chunk start. Absolute event times still reconstruct as `tStart + elapsedEvent`.
- Observer chunk-bridge state is updated once per chunk with `compensatedTimeAdd(...)`; accepted-step elapsed inside the chunk comes from the solver.
- `elapsedbuffer[3]`, `elapsedLastEvent`, and `elapsedThisDown` all remain in the same observation-window-relative coordinate system across continued chunks.
- Trajectory and auxiliary integrals, maxima, minima, and other non-time accumulators remain on the Neumaier or integral helper path.

## Settled public semantics

- `set_tspan((start, end))` resets both the nominal requested window and the hidden solver timebase to `start`.
- `shift_tspan()` advances the nominal requested window by one nominal duration and copies device `tf -> t0`. The next solve uses the attained per-item start time, not `tspan[0]`.
- `shift_x0()` promotes `xf -> x0` and continues the host-side cached problem state.
- `update_x0=True` on `transient()`, `trajectory()`, and `features()` now performs both `shift_x0()` and `shift_tspan()` so repeated solve calls move forward in both state and time by default.
- `get_tspan()` returns the nominal requested window. `get_final_time()` remains the per-item attained absolute final time.

## Evidence

- `transient.cl`, `trajectory.cl`, `features.cl`, and `initializeObserver.cl` now seed each solve from `t0[i]` and rebuild accepted times from the solver-owned chunk-local elapsed pair.
- Fixed-step and adaptive steppers still own chunk-local compensated elapsed and stage-time reconstruction; only the absolute origin moved from shared `tspan[0]` to per-item `t0`.
- `_opencl/buffers.py` and `_opencl/executors.py` now persist device `t0`, reset it on `set_tspan()`, and continue it with `shift_tspan()` via `tf -> t0`.
- `observer_summary.clh` and the semantic event observers now consume solver-owned chunk elapsed directly and persist only chunk-start elapsed offsets with the time-helper family.
- `test/core_numerics/test_features_summary.py::test_rk4_summary_continuation_matches_single_run` still proves split-window continuation for the summary family.
- `test/test_features.py::test_threshold_crossing_continuation_matches_single_run` now proves split-window continuation for a semantic event observer family.
- Focused observer-kernel tests cover threshold, Schmitt, normalized Schmitt, and local-maximum elapsed-buffer behavior after the ownership refactor.
- Seeded stochastic continuation still matches a single run through `test/core_numerics/test_stochastic.py::test_stochastic_euler_ou_continuation_matches_seeded_single_run`.

## Follow-on work

- Internal chunk loops and progress-oriented chunk orchestration should reuse the landed repeated-solve atom instead of adding a second kernel mode.
- Public workflow docs under `docs/` still need a user-facing pass to explain the new default continuation behavior without dragging implementation detail into the docs.
- `get_tspan()` now intentionally reports the nominal requested window rather than per-item attained starts; any future fetchable current-time surface should stay separate from that nominal-window API.

## Public-docs boundary

Public docs should stay usage-focused:

- what continues automatically
- what `set_tspan()` resets
- when to call `shift_tspan()` or `shift_x0()` explicitly
- how `get_tspan()` differs from `get_final_time()`
- what compensated time bookkeeping does and does not fix in float32
