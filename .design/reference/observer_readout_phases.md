# Observer Readout Architecture Phases

Purpose: maintain a short status board for readout-rollout milestones and point to archived rollout history.
Read when: you need to know what parts of the readout rollout are landed versus still open.
Update when: rollout status changes or a previously open readout decision is closed.

## Fast path

- Phases 1-3 (baseline readout rollout) are landed.
- Current active planning is the shared seam for recurring event-trigger controls plus family-local readouts in `.design/next_pr.md`.
- No public family-level readout-selection dataclasses are currently live.
- Historical phase detail and prior implementation checklists are archived in `.design/archived/observer_readout_cleanup_2026_05_31/observer_readout_phases.md`.

## Current status

- Phase 1 (audit and architecture alignment): landed.
- Phase 2 (trajectory stats, period, maxima-count baseline): landed.
- Phase 3 (amplitude rollout across semantic event families): landed.
- Phase 4 (family-alignment decisions): partially open and tracked by the active next-PR target.
- Phase 5 (possible readout-selection surfaces): conditional and not committed.
- Phase 6 (event-state value capture): future/lower priority.

## Remaining open decision set

- Shared seam vs family-local duplication for recurring event-trigger controls plus oscillation-oriented readouts.
- Whether the now-live canonical Schmitt duration outputs (`up duration`, `down duration`, `duty`) should stay family-local or later be folded into any shared oscillation bundle/readout seam.
- Whether any future readout-selection surface is warranted by evidence.

## Where each question lives now

- Active implementation scope and acceptance criteria: `.design/next_pr.md`
- Canonical vs compatibility surface policy: `.design/reference/compatibility_boundary_audit.md`
- Current readout contract and guardrails: `.design/reference/observer_readout_audit.md`
- Deeper observer architecture rationale: `.design/reference/observer_concept_audit.md`
