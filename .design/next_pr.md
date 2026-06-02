# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer readout bundle/selectability proof

## Scope

Resume the threshold-family readout-bundle proof now that solver-owned time-base ownership and continuation semantics are landed. The next narrow target is to define a small set of build-specialized readout bundles on top of the current full-event-geometry defaults without reopening the time-ownership work.

This target should:

- define the first `event geometry`, `oscillation core`, and `trajectory summary` bundles on top of the threshold-family observers
- keep family-specific semantic configs rather than inventing a new cross-family config object
- prove that bundle selection specializes feature schema and persistent layout without regressing existing full-geometry outputs
- keep the landed continuation semantics unchanged and treat the new time-base note as a dependency, not as scope to revisit

**Strictly out of scope:**

- Reopening `set_tspan()` or `shift_tspan()` semantics
- New observer families
- Generic chunk orchestration or progress APIs
- Register-pressure cleanup unrelated to bundle specialization

## Why Now

The continuation blocker is gone. The next highest-value observer work is to make readout scope intentional and testable instead of leaving every event family at one coarse always-on schema.

## Acceptance Criteria

- [ ] Threshold-family observers expose at least the first proved `event geometry`, `oscillation core`, and `trajectory summary` bundles.
- [ ] Bundle choice specializes feature schema and persistent layout through the existing observer-definition and build-key machinery.
- [ ] Full-geometry threshold outputs remain available and backward-compatible.
- [ ] Tests cover schema selection, output naming, persistent-layout specialization, and a no-regression path for the current full threshold outputs.
- [ ] The landed continuation note remains authoritative and no new time-ownership drift appears in docs or code.

## Key Refs

- `.design/reference/observer_readout_audit.md`
- `.design/reference/compatibility_boundary_audit.md`
- `.design/reference/continuation_timebase_note.md`
- `clode/observers/_definitions.py`
- `clode/observers/types.py`
- `clode/simulation/features.py`
- `test/test_features.py`
- `test/test_simulation_contracts.py`
