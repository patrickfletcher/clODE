# Observer Readout Architecture Phases

Purpose: maintain a short status board for readout-rollout milestones and point to archived rollout history.
Read when: you need to know what parts of the readout rollout are landed versus still open.
Update when: rollout status changes or a previously open readout decision is closed.

## Fast path

- Phases 1-3 (baseline readout rollout) are landed.
- Current active planning is the first build-specialized readout-bundle/selectability proof in `.design/next_pr.md`.
- No public family-level readout-selection dataclasses are currently live.
- Historical phase detail and prior implementation checklists are archived in `.design/archived/observer_readout_cleanup_2026_05_31/observer_readout_phases.md`.

## Current status

- Phase 1 (audit and architecture alignment): landed.
- Phase 2 (trajectory stats, period, maxima-count baseline): landed.
- Phase 3 (amplitude rollout across semantic event families): landed.
- Phase 4 (family-alignment decisions): landed far enough to keep family-specific semantic configs and stop chasing a control-only oscillation bundle.
- Phase 5 (possible readout-selection surfaces): now active as a narrow build-specialized proof target.
- Phase 6 (event-state value capture): future/lower priority.

## Remaining open decision set

- Which shared bundles (`event geometry`, `oscillation core`, `trajectory summary`) deserve first build-specialized selection support.
- How Schmitt phase-state extras and neighborhood/local-max family-local outputs should layer on top of those shared bundles.
- Whether any future user-facing readout-selection surface is warranted once a build-specialized proof exists.

## Where each question lives now

- Active implementation scope and acceptance criteria: `.design/next_pr.md`
- Canonical vs compatibility surface policy: `.design/reference/compatibility_boundary_audit.md`
- Current readout contract and guardrails: `.design/reference/observer_readout_audit.md`
- Deeper observer architecture rationale: `.design/reference/observer_concept_audit.md`
