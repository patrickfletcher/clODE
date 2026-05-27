# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer compatibility boundary and bundle decision

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The summary-observer proof slice is landed: `Observer.summary` plus `SummaryObserverSelection` resolve to build-specialized summary variants, while `basic` and `basicall` remain compatibility presets over that same family.
- The threshold catalog now spans four semantic workflows plus one retained legacy full path: `Observer.threshold_crossing`, `Observer.normalized_threshold_crossing`, `Observer.schmitt_trigger`, `Observer.normalized_schmitt_trigger`, and the legacy `Observer.threshold_2` observer.
- `ThresholdCrossingConfig` is shared by the two threshold-crossing families, and `SchmittTriggerConfig` is shared by the two Schmitt families.
- `ObserverParams`, legacy `observer_*` constructor keywords, the retained legacy `threshold_2` observer, and flat barrels remain compatibility surfaces.
- `.design/reference/compatibility_boundary_audit.md` is now the authoritative inventory for which observer-related surfaces are canonical, compatibility-only, or internal-only.

## Why this should be next

The threshold-family work answered the main semantic questions: clODE can keep threshold geometry in runtime settings when the readout schema stays fixed, and it can expose absolute versus warmup-derived parameterization without overloading threshold crossing and Schmitt triggering into one concept.

The remaining pressure point is now the compatibility boundary. Shared config classes mean the selected observer, not the config type alone, determines the semantics, while the older `ObserverParams`, `observer_*`, and alias-heavy surface still remains live. That boundary should be explicit before clODE adds another public reuse layer such as an oscillation-oriented bundle seam.

## Scope

- keep the landed four-family threshold catalog stable while auditing which observer-related surfaces are canonical, compatibility-only, and internal-only
- keep public docs and design notes canonical on semantic names and shared config classes, with compatibility aliases documented only where useful
- use the documented compatibility boundary to decide whether a shared oscillation-oriented bundle seam is justified now or whether family-local duplication should stand for longer
- preserve `ObserverParams`, the current `observer_*` constructor keywords, and flat barrels as compatibility layers rather than expanding them further in this PR

## Key Questions

- Which observer surfaces are canonical, which are compatibility-only, and which should stay internal-only?
- Is the next reusable seam really an oscillation-oriented bundle, or is the current reuse still too thin to justify another public abstraction?
- Should `ObserverParams` and the `observer_*` path remain lightly documented compatibility helpers, or do they need a clearer future support statement?
- How should flat compatibility barrels factor into later public API cleanup once the observer vocabulary settles?

## Design Constraints

- no behavior-breaking public cleanup in this PR
- no broad public config redesign in this PR
- preserve the shared-config rule that the selected observer, not the config class alone, determines absolute versus warmup-derived interpretation
- keep root `.design` docs short and let the detailed alias inventory live in `.design/reference/compatibility_boundary_audit.md`
- keep semantic ownership in Python first and treat `_opencl` as the execution consumer rather than the authoritative definition layer
- do not reopen the threshold naming debate unless the compatibility audit exposes a concrete mismatch that the current semantic names cannot carry

## Non-goals

- no generalized trigger-function or hyperplane DSL in this PR
- no further threshold-family additions unless the audit exposes a concrete missing workflow
- no deprecation campaign yet for compatibility aliases, bundles, or barrels
- no trajectory, solver, or multi-device redesign in this PR

## Suggested Work Slices

1. Maintain the compatibility-boundary reference note and keep docs aligned with it.
2. Decide whether shared oscillation-oriented controls and readouts show enough real reuse to justify a small bundle seam.
3. If not, keep family-local duplication and move the next observer-design effort to a stronger adjacent target.
4. If yes, land the smallest Python-side proof that does not change the executor contract.

## Code-Facing Checklist

- `clode/observers/types.py`, `clode/simulation/features.py`: keep semantic configs canonical and compatibility adapters explicit
- `clode/__init__.py`, `clode/observers/__init__.py`, and the flat compatibility barrels: treat exports as compatibility surfaces unless promoted intentionally
- `.design/reference/compatibility_boundary_audit.md`: authoritative compatibility inventory for observer-related surfaces
- `.design/reference/observer_concept_audit.md`, `docs/feature_extraction.md`, `docs/observers.md`: keep public semantics aligned with the four-family threshold catalog

## Acceptance Criteria

- one authoritative note tracks the current observer compatibility boundary
- user docs prefer semantic names and shared config classes, and mention only the live compatibility aliases
- root planning docs describe the next observer decision as a compatibility-boundary-informed bundle decision rather than another threshold-family addition
- the next implementation candidate is named with evidence rather than inertia

## Follow-on If This Lands Cleanly

1. decide whether one shared oscillation-oriented bundle seam is justified
2. revisit the extremum-family and `nhood1` questions with the compatibility boundary documented
3. defer broader public API cleanup until that next seam is either proven or rejected
