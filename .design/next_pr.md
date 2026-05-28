# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer solution-buffer and state/output bundle audit

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The summary-observer proof slice is landed: `Observer.summary` plus `SummaryObserverSelection` resolve to build-specialized summary variants, while `basic` and `basicall` remain compatibility presets over that same family.
- The semantic trigger catalog now includes `Observer.threshold_crossing`, `Observer.normalized_threshold_crossing`, `Observer.schmitt_trigger`, `Observer.normalized_schmitt_trigger`, `Observer.local_extremum`, and `Observer.neighborhood_return`, while `Observer.threshold_2`, `Observer.local_max`, and `Observer.neighbourhood_2` remain supported retained legacy full workflows.
- The event-method cleanup is now landed for the lean Schmitt families and the retained `threshold_2` timing path: `updateObserverState(...)` advances accepted-step history and continuous reducers, while `eventFunction(...)` plus `computeEventFeatures(...)` own transition detection and refined transition storage.
- The neighborhood-return families already rely on family-local accepted-step buffers to support sampled-anchor detection plus refined exit timing.
- `ObserverParams`, legacy `observer_*` constructor keywords, the retained legacy heavy observers, and flat barrels remain compatibility surfaces, and `.design/reference/compatibility_boundary_audit.md` is the authoritative inventory for that boundary.
- `clode/kernels/clODE_utilities.cl` currently looks more like a private-use helper layer than the main abstraction point for the next observer design decision.

## Why this should be next

The recent event-detection and refinement cleanup resolved the clearest method-role mismatch in the threshold and Schmitt slice. The next source of design friction is not helper placement by itself and not yet the oscillation-oriented readout seam. It is that the live observers still carry several family-local accepted-step buffers and related observer-private state layouts, which makes later observer-state and observer-output bundle work harder to reason about than it needs to be.

A utilities-only PR would mostly reorganize private helpers without answering that structural question. The more useful next decision is whether a shared accepted-step `K`-sample solution-buffer concept is actually the missing prerequisite for clearer observer bundle work, or whether the current family-local buffers are already the right level of duplication.

## Scope

- inventory accepted-step history layouts, buffer-update patterns, and event-local scratch across the observer families that currently depend on recent-step geometry
- decide whether one shared accepted-step `K`-sample solution-buffer concept should precede more observer bundle work
- identify which helper or update patterns should stay family-local even if a shared buffer concept lands
- state whether the oscillation-oriented readout seam should follow immediately after this audit or can proceed without a shared buffer abstraction

## Key Questions

- Which observer families really share the same accepted-step geometry, and which only look similar from a distance?
- Would a shared `K`-sample buffer concept clarify later observer-state and observer-output bundles enough to justify itself now?
- What is the smallest useful abstraction: a conceptual layout contract, shared update helpers, or a concrete kernel-side struct and helper surface?
- Which utilities are still better treated as private caller-side helpers even after that decision?
- If no shared buffer concept lands next, is the oscillation-oriented readout seam still the best immediate follow-on?

## Design Constraints

- no behavior-breaking public cleanup in this PR
- no broad helper-extraction PR just for `clODE_utilities.cl` organization
- keep utilities primarily private-use unless the solution-buffer audit shows a clear shared need
- keep root `.design` docs short and let the detailed alias inventory live in `.design/reference/compatibility_boundary_audit.md`
- keep semantic ownership in Python and the design notes first; kernels consume the decision rather than invent it
- do not widen shared runtime settings unless the audit reveals a real owner boundary that needs them

## Non-goals

- no oscillation-oriented readout seam implementation in this PR
- no generalized trigger-function or hyperplane DSL in this PR
- no deprecation campaign yet for compatibility aliases, bundles, or barrels
- no public fetchable observer-state API yet
- no trajectory, solver, or multi-device redesign in this PR

## Suggested Work Slices

1. Map the accepted-step buffers and event-local scratch each live observer family actually carries.
2. Decide whether one shared solution-buffer concept clarifies later bundle work or would just hide family differences.
3. If yes, name the smallest follow-on proof target and the families it should cover first.
4. If no, document why and re-promote the oscillation-oriented readout seam as the next implementation PR.

## Code-Facing Checklist

- `clode/kernels/observers/*.clh`: inventory accepted-step history, event-local scratch, and duplicated buffer updates without widening behavior in the audit PR
- `clode/kernels/clODE_utilities.cl`: treat helper placement as secondary to the solution-buffer decision; avoid turning private-use helpers into a premature public kernel abstraction
- `.design/reference/observer_concept_audit.md`, `.design/reference/semantic_layout_audit.md`, `.design/reference/compatibility_boundary_audit.md`: keep the trigger-story, state/bundle question, and compatibility boundary aligned
- `.design/package_state.md`, `.design/ideas.md`: keep the active target and dependency ordering explicit

## Acceptance Criteria

- one authoritative design note states whether a shared accepted-step solution-buffer concept is a prerequisite for clearer observer bundle work
- root planning docs agree that utilities alone are not the next proof target
- the dependency between solution-buffer work and oscillation-oriented bundles is explicit in `ideas.md` and `next_pr.md`
- the first follow-on implementation target after the audit is named clearly

## Follow-on If This Lands Cleanly

1. If a shared buffer concept is justified, land the smallest solution-buffer or shared-update-helper proof before resuming observer bundle design.
2. If it is not justified, return to the oscillation-oriented readout seam with the duplication rationale documented.
3. After that, measure heavy observer-state footprint and register pressure against the chosen direction.
