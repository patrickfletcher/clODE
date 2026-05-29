# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer docs and planning cleanup

## Assumed Repo State

- The observer consolidation work is landed in code and covered by focused tests.
- `Observer.local_max` is canonical; `Observer.local_extremum` remains a compatibility alias.
- Threshold families always expose event times and event count.
- Event stream naming is standardized in the live observer surface.

## Why this should be next

The live observer surface is now stable enough that the remaining work is to make the root docs and planning board tell the same story without pilot-era wording or duplicate stale notes.

## Scope

- Remove stale docs language that mentions local-extremum polarity selection.
- Remove stale docs/design language that mentions `ThresholdCrossingReadoutSelection`.
- Ensure event naming examples in docs match live feature names.
- Trim `ideas.md` and related root planning notes so they stay open-item only and do not repeat landed observer history.

## Non-goals

- No new observer families.
- No new readout-selection APIs.
- No semantic changes to trigger geometry.

## Code-Facing Checklist

- `docs/feature_extraction.md`: local_max canonical wording, no polarity or readout-selection pilot text.
- `docs/observers.md`: readout inventory and event-stream naming align with live outputs.
- `.design/package_state.md`: snapshot text reflects removed pilot surface.
- `.design/ideas.md`: remove or archive stale local-extremum-decision line if no longer open.

## Acceptance Criteria

- Docs no longer reference `ThresholdCrossingReadoutSelection` or `ExtremumPolarity`.
- Docs describe `Observer.local_max` as canonical; `local_extremum` appears as compatibility-only where needed.
- Event timestamp naming examples match live feature names.
- Root `.design` docs stay terse, current, and free of completed observer-implementation history.
- Focused observer contract tests pass.

## Follow-on If This Lands Cleanly

1. Add optional period/amplitude summary outputs to `Observer.summary`.
2. Revisit full observer config/readout seams only after standardized observer set is stable.
