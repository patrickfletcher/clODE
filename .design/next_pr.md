# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Oscillation-oriented observer bundle seam: decide and first implementation slice

## Scope

Determine how recurring event-trigger-family controls (`min_amp`, `max_event_count`) and observer-family readouts (for example Schmitt `up duration`, `down duration`, and `duty`) should be packaged through a shared seam versus family-local composition across event-observer families.

Resolved policy inputs for this PR target:

- Schmitt duration readouts (`up duration`, `down duration`, `duty`) are core canonical Schmitt measurements.
- Semantic config surfaces should expose `feature_var` whenever kernels consume both `eVarIx` and `fVarIx` semantics.
- `max_event_count` is a canonical event-trigger-family limiter and early-stop control rather than merely a legacy compatibility field.
- `min_amp` is a user-specified event gate on variation in `event_var`; current family differences are about how that variation is measured, not about its purpose.
- The current event-triggering semantic configs now all expose `min_amp` and `max_event_count`; the remaining question is packaging, not whether those controls belong semantically.
- `eps_dx` currently has no semantic-family use case and should remain compatibility-only unless a later audit finds one.

Land the decision as either a first shared seam or an explicit kept-separate record with rationale, and update the relevant docs and planning notes to reflect it.

**Strictly out of scope:**

- Observer-state or register-pressure cleanup (later pass)
- New observer families beyond what the bundle decision requires
- Compatibility-barrel retirement
- Numerical evidence or publication work

## Why Now

The solution-buffer audit is settled, K=2/K=3 shared helpers are live across canonical event-observer families, and the `compatibility_boundary_audit.md` already surfaces the bundle question as an open follow-up. The decision about whether to share or separate the oscillation-oriented seam is what determines the authoring and docs story for all future observer additions.

## Acceptance Criteria

- [x] `min_amp` and `max_event_count` are surveyed across the current event-triggering families, their live semantics and current adoption scope are documented, and at least one candidate additional family-local readout is surveyed alongside them.
- [x] Canonical Schmitt readout direction is documented explicitly: `up duration`, `down duration`, and `duty` are core Schmitt outputs and are now exposed on canonical Schmitt schemas.
- [x] Config-surface direction is documented explicitly: threshold, Schmitt, and neighborhood families that split trigger geometry from extrema/amplitude behavior expose `feature_var`.
- [ ] A decision is recorded in `.design/reference/compatibility_boundary_audit.md`: either (a) a shared oscillation-bundle seam design with a named config object and at least one family adoption, or (b) an explicit keep-separate rationale explaining why family-local duplication is preferable.
- [ ] If a shared seam is chosen: at least one event-observer family adopts it, a test covers the new config surface, and docs reflect the new preferred authoring pattern.
- [ ] If keep-separate is chosen: `ideas.md` oscillation-bundle item is updated to reflect the decision, and the `compatibility_boundary_audit.md` bundle question is closed with rationale.
- [ ] `package_state.md` is updated to reflect whichever direction lands.
- [ ] No regressions in `test/test_features.py`, `test/test_simulation_contracts.py`, or `test/kernel_components/`.

## Key Refs

- `.design/reference/observer_concept_audit.md` — deeper observer rationale and bundle design axes
- `.design/reference/compatibility_boundary_audit.md` — current canonical/compatibility boundary and open bundle question
- `.design/tmp/observer_threshold_2.clh` — Schmitt duration/duty prototype reference for canonical-surface integration work
- `clode/observers/types.py`, `clode/simulation/features.py` — Python observer surface
- `docs/feature_extraction.md`, `docs/observers.md` — user-facing observer docs
