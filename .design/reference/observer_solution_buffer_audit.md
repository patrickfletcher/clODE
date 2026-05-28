# Observer Solution-Buffer Audit

Purpose: decide whether a shared accepted-step solution-buffer concept should precede observer state/output bundle work, and record the minimum contract that current observer families actually share.
Read when: planning observer-state or observer-output bundle work, deciding whether to centralize accepted-step history updates, or scoping follow-on observer kernel cleanup.
Update when: the shared-buffer decision changes, a follow-on proof lands, or observer families materially change their accepted-step history layouts.

## Bottom line

- A shared accepted-step `K`-sample solution-buffer concept should come before further observer bundle work.
- The shared concept should start as a layout-and-update contract, not a single forced `ObserverState` struct.
- The strongest overlap is in accepted-step history advancement (`t`, `x`, and sometimes `dx`) for `K=2` and `K=3` families.
- Event refinement guidance should default to linear or inverse-linear threshold timing where already proven stable, keep bounded three-sample quadratic extrema for current extremum families, and treat cubic Hermite refinement as a targeted follow-on where endpoint slopes are already reliable.
- Heavy legacy observers still need family-local state and bundle logic beyond that shared history core.
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

### Families that already match a common accepted-step history core

- `observer_threshold_crossing.clh`: `K=2` with `tbuffer[2]` and one event-variable sample buffer `xbuffer[2]`.
- `observer_normalized_threshold_crossing.clh`: same `K=2` core as absolute threshold crossing.
- `observer_schmitt_trigger.clh`: `K=2` with `tbuffer[2]` and event-variable `xbuffer[2]` plus Schmitt state.
- `observer_normalized_schmitt_trigger.clh`: same `K=2` core with warmup-derived thresholds.
- `observer_neighborhood_return.clh`: `K=2` with `tbuffer[2]` and full-state `xbuffer[2 * N_VAR]` for normalized distance checks.
- `observer_local_extremum.clh`: `K=3` with `tbuffer[3]`, `xbuffer[3]`, and `dxbuffer[3]` for three-sample extremum refinement.

### Families that share history mechanics but carry larger legacy bundles

- `observer_threshold_2.clh`: `K=3` over full state and slope (`xbuffer[3 * N_VAR]`, `dxbuffer[3 * N_VAR]`) plus elapsed-time buffers, periodic summaries, and wide trajectory/aux outputs.
- `observer_local_maximum.clh`: `K=3` over full state and slope plus elapsed-time buffers and heavier summary/event lists.
- `observer_neighborhood_2.clh`: `K=3` over full state and slope plus warmup ranges, center-point state, and broader period/summary bundles.
- `observer_neighborhood_1.clh`: `K=3` over full state and slope with mutable normalization and broader legacy summary state.

### Families where a `K`-sample event-history concept is not the main abstraction

- `observer_summary.clh`, `observer_basic.clh`, `observer_basic_allVar.clh`: online reducers with no event detector history window; they should remain outside the shared accepted-step history proof.

## What is shared versus local

### Shared candidate surface

- Accepted-step index advancement for `K=2` and `K=3` time buffers.
- Accepted-step index advancement for aligned state and slope buffers when present.
- Small interpolation-ready local sample views consumed by event or refinement helpers.

### Scope and ownership for the shared buffer contract

- The accepted-step solution buffer should remain observer-owned runtime state (`ObserverState`) so continuation semantics stay intact and family state can be persisted through `observer_states` global storage.
- `features.cl` should stay the orchestration kernel that advances the solver and calls observer methods; it should not own per-family rolling buffer layout policy.
- Shared history shift/update helpers should live with observer-kernel support code (observer headers and possibly narrowly scoped shared helper utilities), with explicit address-space contracts and no hidden generic-pointer assumptions.
- `clODE_utilities.cl` can host math primitives used by refinement routines, but observer-family policy and layout-specific transitions stay in observer code.

### Keep family-local

- Warmup-derived threshold/range learning and detector initialization policy.
- Event semantics (`threshold`, `Schmitt`, `extremum`, `neighborhood`) and state-machine logic.
- Heavy legacy bundle accumulation (period, duty, amplitude, aux/state aggregate outputs).
- Family-specific retained sparse output arrays and continuation-reset policy.

## Interpolation and event-refinement guidance

This section incorporates the recommendations in `ode_event_interpolation_note.md` and binds them to the current observer implementation path.

- Keep detection and refinement distinct: detection remains family semantics, refinement improves stored event time or event state after detection.
- Threshold crossings remain the easiest and best-conditioned refinement target; current linear or inverse-linear crossing time helpers remain the baseline.
- Extremum refinement remains intrinsically harder; retain bounded three-sample quadratic extrema as the stable default for current extremum families.
- Where endpoint slope samples are already present and reliable for the event variable, cubic Hermite refinement is a valid follow-on option for selected families, not a cross-family requirement for slice 2.
- Do not widen interpolation order globally in this proof; higher-order dense output belongs to a separate numerics-driven decision.

## Relation to future user-selectable feature bundles

- A shared accepted-step buffer concept should make future bundle selection cleaner by separating geometry acquisition (history updates) from bundle policy (which summaries or periodic outputs to accumulate).
- Future selectable bundles (for example summary aggregates, periodic features, and duty-cycle groups) should remain a separate declaration and schema decision on the Python and observer-definition side.
- Slice 2 should prepare this seam by reducing duplicate history-update code, not by implementing bundle selection itself.

## Decision

A shared accepted-step solution-buffer concept is a prerequisite for clearer observer bundle work, but only as a narrow core contract.

The contract should state:

- `K` is family-specific (`K=2` or `K=3` in current families).
- The shared layer owns accepted-step history shift/update mechanics.
- Families choose whether the shared history carries event-variable-only samples or full `N_VAR`/`N_AUX` views.
- Families retain full ownership of event semantics, warmup policy, and readout bundles.

## Recommended follow-on proof

1. Add shared accepted-step history update helpers for `K=2` and `K=3` usage patterns in observer code paths.
2. Apply the helpers first to one lean family (`threshold_crossing` or `schmitt_trigger`) and one heavier family (`threshold_2` or `neighborhood_2`) to test both ends.
3. Keep helper scope to history shift/update and interpolation-ready sample exposure only; do not centralize event semantics or bundle accumulation in the same PR.
4. Measure readability and duplication reduction in touched families before broad rollout.

## Non-goals for the proof

- No immediate unification of all observer state structs.
- No broad `clODE_utilities.cl` reorganization unrelated to accepted-step history updates.
- No behavior or public-schema changes for existing observer outputs.
