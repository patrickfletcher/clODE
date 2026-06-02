# Observer Solution-Buffer Audit

Purpose: decide whether a shared accepted-step solution-buffer concept should precede observer state/output bundle work, and record the minimum contract that current observer families actually share.
Read when: planning observer-state or observer-output bundle work, deciding whether to centralize accepted-step history updates, or scoping follow-on observer kernel cleanup.
Update when: the shared-buffer decision changes, a follow-on proof lands, or observer families materially change their accepted-step history layouts.

## Bottom line

- A shared accepted-step `K`-sample solution-buffer concept should come before further observer bundle work.
- The shared concept should start as a layout-and-update contract, not a single forced `ObserverState` struct.
- The strongest overlap is in accepted-step history advancement for solver-relative elapsed time plus aligned `x`, `dx`, and optional `aux` samples across the current `K=3` families.
- The semantic event observers no longer need mirrored absolute-time history buffers; they now retain a start-time origin plus solver-relative elapsed history and reconstruct absolute event timestamps when events are stored.
- Event refinement guidance should default to linear or inverse-linear threshold timing where already proven stable, keep bounded three-sample quadratic extrema for current extremum families, and prioritize three-point quadratic threshold-crossing time/state refinement as the higher-order path.
- Feature-rich observer families still need family-local state and bundle logic beyond that shared history core.
- For the active observer path, the work-item-local accepted-step buffer belongs to observer-owned state and observer helper contracts, not as an ad hoc local buffer in `features.cl`.
- The smallest follow-on proof target is shared accepted-step history update helpers for `K=2` and `K=3`, while retaining family-local event semantics and feature accumulation.

## Fast path

Stop after this section unless you are implementing the follow-on helper proof or reassessing family layouts.

- This note is the authoritative decision for whether shared accepted-step history should precede additional observer bundle work.
- This note also records the current scope decision for where accepted-step buffers should live in kernel code.
- Use `## Family inventory` to see what is truly shared versus still family-local.
- Use `## Interpolation and event-refinement guidance` before widening any event-time or event-state interpolation helper surface.
- Use `## Recommended follow-on proof` for the next implementation target and non-goals.
- Use root `.design` docs for active sequencing and PR scope.

## Family inventory

After the readout rollout, all semantic event families adopted the shared K=3 accepted-step history update helpers. The earlier lean-shell design used event-variable-only K=2 buffers for threshold and Schmitt families; that design no longer reflects the live kernels.

### Families that match the common K=3 accepted-step history core

- `observer_threshold_crossing.clh`: `K=3` with `elapsedbuffer[3]`, `xbuffer[3 * N_VAR]`, `dxbuffer[3 * N_VAR]`, and `auxbuffer[3 * N_AUX]`; event detection reads only the `eVarIx` slice at positions 1–2 while stored event timestamps are rebuilt from `tStart + elapsed`.
- `observer_normalized_threshold_crossing.clh`: same `K=3` elapsed-buffer layout as absolute threshold crossing.
- `observer_schmitt_trigger.clh`: `K=3` with `elapsedbuffer[3]`, `xbuffer[3 * N_VAR]`, `dxbuffer[3 * N_VAR]`, and `auxbuffer[3 * N_AUX]`; Schmitt state transitions read the `eVarIx` slice at positions 1–2 and store full up/down transition geometry.
- `observer_normalized_schmitt_trigger.clh`: same `K=3` elapsed-buffer layout as absolute Schmitt, with warmup-derived thresholds.
- `observer_local_maximum.clh`: `K=3` with `elapsedbuffer[3]`, `xbuffer[3 * N_VAR]`, `dxbuffer[3 * N_VAR]`, and `auxbuffer[3 * N_AUX]` for three-sample extremum detection and refinement.
- `observer_normalized_neighborhood_return.clh`: `K=3` with `elapsedbuffer[3]`, `xbuffer[3 * N_VAR]`, `dxbuffer[3 * N_VAR]`, and `auxbuffer[3 * N_AUX]`; the `elapsedbuffer` tracks intra-period elapsed time for anchor interpolation and full event geometry.

### Families that carry additional per-family state beyond the K=3 history core

The shared K=3 layout is not the full picture for heavier families; these carry per-family readout state that remains local:

- `observer_normalized_schmitt_trigger.clh`: warmup-derived threshold state (`xUp`, `xDown`) plus period, amplitude, and duration bundle fields.
- `observer_local_maximum.clh`: sparse max/min event geometry lists, full trajectory-stat fields (`xTrajectoryMax`, `xTrajectoryMin`, etc.), and period/amplitude bundle.
- `observer_normalized_neighborhood_return.clh`: center-point `x0`, warmup range state, per-event elapsed bookkeeping, and full trajectory-stat fields with period/amplitude bundle.

### Families where a `K`-sample event-history concept is not the main abstraction

- `observer_summary.clh`: online reducers with no event detector history window; they should remain outside the shared accepted-step history proof.

## What is shared versus local

### Shared candidate surface

- Accepted-step index advancement for solver-relative elapsed buffers.
- Accepted-step index advancement for aligned state and slope buffers when present.
- Accepted-step index advancement for aligned auxiliary buffers when a family retains full event geometry.
- Small interpolation-ready local sample views consumed by event or refinement helpers.

### Scope and ownership for the shared buffer contract

- The accepted-step solution buffer should remain observer-owned runtime state (`ObserverState`) so continuation semantics stay intact and family state can be persisted through `observer_states` global storage.
- `features.cl` should stay the orchestration kernel that advances the solver and calls observer methods; it should not own per-family rolling buffer layout policy.
- Shared history shift/update helpers should live with observer-kernel support code (observer headers and possibly narrowly scoped shared helper utilities), with explicit address-space contracts and no hidden generic-pointer assumptions.
- `clODE_utilities.cl` can host math primitives used by refinement routines, but observer-family policy and layout-specific transitions stay in observer code.

### Keep family-local

- Warmup-derived threshold/range learning and detector initialization policy.
- Event semantics (`threshold`, `Schmitt`, `extremum`, `neighborhood`) and state-machine logic.
- Family-specific bundle accumulation (period, duration, duty, amplitude, aux/state aggregate outputs).
- Family-specific retained sparse output arrays and continuation-reset policy.

## Interpolation and event-refinement guidance

This section incorporates the recommendations in `ode_event_interpolation_note.md` and binds them to the current observer implementation path.

- Keep detection and refinement distinct: detection remains family semantics, refinement improves stored event time or event state after detection.
- Threshold crossings remain the easiest and best-conditioned refinement target; current linear or inverse-linear crossing time helpers remain the baseline.
- Extremum refinement remains intrinsically harder; retain bounded three-sample quadratic extrema as the stable default for current extremum families.
- For threshold crossing time/state refinement beyond inverse-linear timing, evaluate three-point quadratic interpolation on accepted-step buffers as the preferred higher-order path.
- Do not widen interpolation order globally in this proof; higher-order dense output belongs to a separate numerics-driven decision.

## Relation to future user-selectable feature bundles

- A shared accepted-step buffer concept should make future bundle selection cleaner by separating geometry acquisition (history updates) from bundle policy (which summaries or periodic outputs to accumulate).
- Future selectable bundles (for example summary aggregates, periodic features, and duty-cycle groups) should remain a separate declaration and schema decision on the Python and observer-definition side.
- Slice 2 should prepare this seam by reducing duplicate history-update code, not by implementing bundle selection itself.

## Decision

A shared accepted-step solution-buffer concept is a prerequisite for clearer observer bundle work, but only as a narrow core contract.

**Status**: The contract is established and the follow-on proof is landed.

The K=2 and K=3 shared accepted-step history update helpers (`advanceAcceptedStepHistory2`, `advanceAcceptedStepHistory2ByVariable`, `advanceAcceptedStepHistory3`, `advanceAcceptedStepHistory3ByVariable`) are now live in `clode/kernels/observers.cl` and adopted by the semantic observer families: `threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, `local_max`, and `normalized_neighborhood_return`. Family event semantics and bundle accumulation remain local.

**The contract states**:

- `K` is family-specific (`K=2` or `K=3` in current families).
- The shared layer owns accepted-step history shift/update mechanics.
- Families choose whether the shared history carries event-variable-only samples or full `N_VAR`/`N_AUX` views.
- Families retain full ownership of event semantics, warmup policy, and readout bundles.

**Next planning target**: a build-specialized readout bundle/selectability proof on top of the now-standardized history and compensated-mean helpers. See `next_pr.md` and `.design/ideas.md`.
