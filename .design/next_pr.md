# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Oscillation-oriented observer readout seam

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The summary-observer proof slice is landed: `Observer.summary` plus `SummaryObserverSelection` resolve to build-specialized summary variants, while `basic` and `basicall` remain compatibility presets over that same family.
- The semantic trigger catalog now includes `Observer.threshold_crossing`, `Observer.normalized_threshold_crossing`, `Observer.schmitt_trigger`, `Observer.normalized_schmitt_trigger`, `Observer.local_extremum`, and `Observer.neighborhood_return`, while `Observer.threshold_2`, `Observer.local_max`, and `Observer.neighbourhood_2` remain supported retained legacy full workflows.
- `ThresholdCrossingConfig`, `SchmittTriggerConfig`, `LocalExtremumConfig`, and `NeighborhoodReturnConfig` are live semantic config surfaces; the selected observer still owns warmup-derived versus absolute interpretation where a config class is shared.
- `ObserverParams`, legacy `observer_*` constructor keywords, the retained legacy heavy observers, and flat barrels remain compatibility surfaces.
- `.design/reference/compatibility_boundary_audit.md` is now the authoritative inventory for which observer-related surfaces are canonical, compatibility-only, or internal-only.

## Why this should be next

The trigger-semantics work answered the main naming and family-boundary questions: clODE can keep threshold, extremum polarity, and neighborhood-return trigger geometry separate while still exposing lean semantic observers next to retained heavier legacy workflows.

The remaining observer-design pressure point is no longer naming. It is that the heavier legacy observers still bundle trigger semantics, oscillation readouts, and persistent-state size together. The next decision should therefore be whether recurring oscillation-oriented readouts and controls deserve a small shared declaration seam, or whether they should stay family-local until there is stronger evidence.

## Scope

- audit recurring oscillation-oriented readouts and control knobs across `local_max`, `threshold_2`, and `neighbourhood_2`
- decide which of those outputs or controls justify a shared selectable readout seam and which should stay family-local
- keep the retained heavy observers supported while clarifying that the lean semantic observers are the preferred default workflows
- identify the smallest Python-side proof target for state-size control that does not widen shared runtime settings prematurely

## Key Questions

- Which readouts recur strongly enough across the retained heavy observers to justify one shared seam?
- Should selected readouts change build-specialized layout, runtime policy, or both?
- How should `max_event_timestamps` and selected readouts combine to control observer-state size without obscuring the current semantic observers?
- What remains compatibility-only after that seam, and what becomes canonical?

## Design Constraints

- no behavior-breaking public cleanup in this PR
- no trigger-semantics renaming pass in this PR
- preserve the shared-config rule that the selected observer, not the config class alone, determines absolute versus warmup-derived interpretation
- keep root `.design` docs short and let the detailed alias inventory live in `.design/reference/compatibility_boundary_audit.md`
- keep semantic ownership in Python first and treat `_opencl` as the execution consumer rather than the authoritative definition layer
- do not widen shared runtime settings unless the selected readout seam really needs it

## Non-goals

- no generalized trigger-function or hyperplane DSL in this PR
- no deprecation campaign yet for compatibility aliases, bundles, or barrels
- no public fetchable observer-state API yet
- no trajectory, solver, or multi-device redesign in this PR

## Suggested Work Slices

1. Inventory recurring heavy-observer readouts and control knobs.
2. Decide whether one oscillation-oriented readout seam is justified now.
3. If yes, land the smallest Python-side declaration proof with matching docs and tests.
4. If not, document why family-local duplication should stand and shift the next effort to state-footprint measurement instead.

## Code-Facing Checklist

- `clode/observers/types.py`, `clode/observers/_definitions.py`, `clode/simulation/features.py`: keep semantic configs canonical and any new readout-selection seam Python-owned
- `clode/kernels/observers/*.clh`: treat kernel state size as a consequence of the selected declaration rather than the vocabulary source
- `.design/reference/compatibility_boundary_audit.md`, `.design/reference/observer_concept_audit.md`: keep the semantic-versus-legacy observer story and the next readout question aligned
- `docs/feature_extraction.md`, `docs/observers.md`: keep user docs focused on semantic observer workflows first and heavier legacy bundles second

## Acceptance Criteria

- one authoritative design note states whether a shared oscillation-oriented readout seam is justified
- the active design language distinguishes trigger semantics from optional heavier readout bundles
- any selected proof target names how observer-state size will shrink or stay bounded
- docs and root planning notes agree on the next observer decision

## Follow-on If This Lands Cleanly

1. measure observer-state footprint and register pressure for the retained heavy kernels against the selected readout model
2. decide whether `neighbourhood_1` still deserves to stay public once the heavier readout story is clearer
3. defer broader public API cleanup until the readout seam is either proven or rejected
