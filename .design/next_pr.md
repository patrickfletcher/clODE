# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Warmup-derived fractional threshold crossing

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The summary-observer proof slice is now landed: `Observer.summary` plus `SummaryObserverSelection` resolve to build-specialized summary variants, while `basic` and `basicall` remain compatibility presets over that same family.
- The simulation state and output ownership pass is wrapped tightly enough that integration settings, trajectory output policy, observer runtime settings, event-output policy, solver state, persistent observer state, and fetched outputs have clearer named homes on the Python side.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.
- The first deterministic static-trigger slice is now landed: `Observer.threshold_crossing` is the preferred semantic alias for the current one-pass absolute threshold observer, while `threshold_1` remains the compatibility spelling.
- `Observer.schmitt_trigger` is now the preferred semantic alias for the current warmup-derived fractional Schmitt family, while `threshold_2` remains the compatibility spelling.
- The threshold-family semantic config seam is now also landed for the current built-ins: `ThresholdCrossingConfig` and `SchmittTriggerConfig` expose family-relevant knobs while the broad compatibility bundle remains intact.
- The next decision is therefore narrower again: use that landed seam to add the missing warmup-derived fractional directional-threshold family.

## Why this should be next

The summary slice answered the first implementation question, and the landed `threshold_1` slice answered the next counterexample: a useful event family can keep threshold value and crossing direction as runtime settings when its readout schema stays fixed.

The follow-on design evidence is now stronger than a generic “threshold versus hysteresis” framing. `threshold_2` already shows that warmup-derived fractional thresholds are useful when event-variable amplitudes shift across a sweep, while `threshold_1` shows that a lean directional threshold family can stay runtime-only.

The immediate architectural bottleneck is no longer the threshold-family config seam itself. That seam now exists on the Python side, and the executor still only consumes normalized runtime settings and event-output settings.

The right next step is therefore back to the missing family that motivated the re-audit: a lean warmup-derived fractional directional-threshold observer with threshold-style timestamps and count, separate from the heavier Schmitt-style outputs.

## Scope

- keep simple directional threshold crossing and Schmitt triggering as separate public concepts even though one is a degenerate implementation case of the other
- use the landed threshold-family config seam instead of extending the generic threshold compatibility bundle again
- add the missing warmup-derived directional-threshold family with a semantic config that matches threshold crossing rather than Schmitt triggering
- keep `threshold_2` documented as the warmup-derived Schmitt-style family and bound the role of slope gates to noise-oriented workflows
- preserve `ObserverParams` and the current `observer_*` constructor keywords as compatibility layers rather than expanding them further
- keep `nhood2` as a specialized periodicity detector rather than a template and keep `nhood1` under explicit keep-or-retire review
- keep generalized scalar trigger functions or hyperplane crossings as a later generalization rather than the immediate implementation surface
- keep trajectory variable-subset work and any public custom-observer DSL out of this PR

## Key Questions

- Which semantic name should the missing warmup-derived directional-threshold family use so it stays distinct from `Observer.threshold_crossing` and `Observer.schmitt_trigger`?
- Which readouts belong to that lean fractional directional-threshold family, and which belong only to the heavier Schmitt-style family?
- Should the new family keep only one threshold scalar in its semantic config object, or should it reserve room for later threshold-conditioned readouts on `fVarIx`?
- How long should `fVarIx` remain only a preserved role in the lean threshold family before a richer event-conditioned readout is justified?
- Should `localmax` and mirrored `localmin` be described as one future extremum family with a polarity selector?
- Is `nhood1` worth keeping once the first deterministic static-trigger slice is landed, or should it move toward retirement unless a stronger workflow is documented?

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, solver-owned status boundary, continuation semantics, and current proof layer
- do not reintroduce solver-owned step, status, or time diagnostics through observer outputs or observer-private state
- keep semantic ownership in Python first and treat `_opencl` as the execution consumer rather than the authoritative definition layer
- no user-facing custom-observer DSL, broad inheritance hierarchy, or full code-conversion surface in the same PR
- avoid turning this PR into a benchmark, register-pressure, or "add every missing observer" campaign; the point is to record and prove the right deterministic event-family boundary before broadening the catalog
- preserve the landed summary-family behavior and current public presets while evaluating follow-on directions
- keep public docs and design notes explicit about what is clarified here and what remains deferred
- preserve the landed lean `threshold_1` schema instead of expanding it opportunistically without a fresh design pass
- do not collapse `threshold_1` and `threshold_2` into one overloaded public option just because the underlying trigger logic can share machinery
- do not let the current compatibility bundle dictate the semantic knob names for future observer-family config objects
- do not reopen the generic threshold naming debate inside this PR unless the new family truly forces a better explicit qualifier than the current semantic vocabulary supports

## Non-goals

- no implicit or IMEX solver work
- no multi-device work
- no broader diverged-time continuation-policy redesign or matched device-side current-time model in the same PR
- no public bundled solver-stats object in the same PR
- no solver-family or stepper-model redesign in the same PR
- no public observer-state fetch API or public rename of `ObserverOutput` in the same PR
- no trajectory variable-subset implementation in the same PR
- no generalized trigger-function or hyperplane DSL in the same PR
- no citation metadata, release-tag, or broader repo-surface cleanup in the same PR
- no broad public-config rewrite beyond the threshold-family seam needed for this slice
- no broader observer-family config rollout beyond the threshold catalog in the same PR

## Suggested Work Slices

1. Define the semantic config and readout shape for the lean warmup-derived directional-threshold family.
2. Land that family through the existing threshold-family config seam instead of extending `ObserverParams` again.
3. Update docs and tests so the new family is clearly distinct from both `Observer.threshold_crossing` and `Observer.schmitt_trigger`.
4. If one ambiguity remains after that, use a very small trial artifact instead of widening production code immediately.

## Code-Facing Checklist

- `clode/observers/_definitions.py`, `clode/observers/types.py`, `clode/observers/metadata.py`: preserve the separate trigger-variable and feature-variable roles while shaping the next event-family declaration layer
- `clode/_opencl/source_builder.py`, `clode/_opencl/models.py`, `clode/_opencl/executors.py`: document how build signature, source preamble, and runtime-setting uploads currently divide responsibility
- `clode/observers/types.py`, `clode/simulation/features.py`: keep the landed threshold-family config seam stable while adding the next threshold family through it
- `clode/kernels/observers/observer_threshold_1.clh`, `clode/kernels/observers/observer_threshold_2.clh`: remain the semantic evidence base for the missing fractional directional-threshold family
- `.design/reference/observer_concept_audit.md`, `.design/ideas.md`, `.design/development_roadmap.md`, `.design/package_state.md`: keep the planning and current-state narrative aligned with the landed config seam and the next missing threshold family

## Acceptance Criteria

- the design docs and public docs distinguish threshold crossing from Schmitt triggering as separate user-facing concepts
- the docs make clear that `threshold_2` currently uses warmup-derived fractional thresholds and Schmitt-style hysteresis, with slope gates as a narrower noise-oriented refinement
- the threshold-family semantic config seam remains the preferred public surface for threshold observers, while the compatibility bundle continues to work
- the next implementation candidate is named with evidence: a lean two-pass fractional directional-threshold family that lands through the semantic threshold config seam
- the design docs make clear which threshold-family follow-ons still fit the lean fixed schema and which knobs belong only to the Schmitt-style family
- the docs record whether mirrored `localmin` belongs in a future extremum family and whether `nhood1` is a keep-or-retire case
- generalized scalar trigger functions or hyperplane crossings remain explicitly deferred rather than drifting back into near-term scope

## Follow-on If This Lands Cleanly

1. implement the lean warmup-derived fractional directional-threshold family through the landed semantic threshold config seam
2. decide when threshold-conditioned readouts become rich enough to require declaration-level schema control, and which of those belong only to Schmitt-style families
3. revisit generalized hyperplane triggers only after the landed scalar-crossing family and its next adjacent extension show what reusable trigger packaging is still missing
