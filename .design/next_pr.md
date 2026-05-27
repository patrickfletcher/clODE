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
- The first deterministic static-trigger slice is now landed: `Observer.threshold_1` provides one-pass absolute threshold crossings with runtime direction selection and a lean threshold-event timestamp stream.
- `Observer.threshold_2` remains valuable because it interprets its threshold inputs as warmup-derived fractions of the observed event-variable amplitude and then applies Schmitt-style up/down state with optional slope gates.
- The next decision is therefore narrower than before: separate threshold topology from threshold parameterization cleanly and fill the missing warmup-derived directional-threshold corner without losing the lean boundary `threshold_1` just proved.

## Why this should be next

The summary slice answered the first implementation question, and the landed `threshold_1` slice answered the next counterexample: a useful event family can keep threshold value and crossing direction as runtime settings when its readout schema stays fixed.

The follow-on design evidence is now stronger than a generic “threshold versus hysteresis” framing. `threshold_2` already shows that warmup-derived fractional thresholds are useful when event-variable amplitudes shift across a sweep, while `threshold_1` shows that a lean directional threshold family can stay runtime-only.

The missing proof target is therefore a lean two-pass directional-threshold family whose threshold is expressed as a warmup-derived fraction instead of an absolute state value. That target cleanly separates threshold parameterization from Schmitt topology. It is a more informative next slice than adding hysteresis to `threshold_1`, and it leaves the current `threshold_2` free to stay the Schmitt-style family rather than the catch-all threshold bucket.

## Scope

- keep simple directional threshold crossing and Schmitt triggering as separate public concepts even though one is a degenerate implementation case of the other
- make the threshold taxonomy explicit: absolute versus warmup-derived parameterization, and single-boundary crossing versus Schmitt-style hysteresis
- shape the next implementation slice around a lean two-pass fractional directional-threshold family with threshold-style readouts rather than period-style Schmitt outputs
- keep `threshold_2` documented as the warmup-derived Schmitt-style family and bound the role of slope gates to noise-oriented workflows
- keep `nhood2` as a specialized periodicity detector rather than a template and keep `nhood1` under explicit keep-or-retire review
- keep generalized scalar trigger functions or hyperplane crossings as a later generalization rather than the immediate implementation surface
- keep trajectory variable-subset work and any public custom-observer DSL out of this PR

## Key Questions

- Should the public observer catalog keep separate threshold and Schmitt trigger options even though they overlap mathematically?
- How should warmup-derived fractional thresholds be exposed without making simple threshold crossing inherit Schmitt-only concepts?
- Which readouts belong to a lean fractional directional-threshold family, and which belong only to the heavier Schmitt-style family?
- Should slope thresholds remain tied to noisy or stochastic workflows rather than appear as baseline threshold knobs?
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
- no immediate public rename of `threshold_2` in the same PR; naming cleanup can wait until the threshold taxonomy is better represented in the built-in catalog

## Suggested Work Slices

1. Update the observer audit and user docs so they distinguish absolute threshold crossing, warmup-derived fractional thresholds, and Schmitt-style hysteresis explicitly.
2. Define the lean fractional directional-threshold family: one warmup-derived threshold, direction selection, preserved `eVarIx` versus `fVarIx`, and a threshold-style timestamp-plus-count readout.
3. Keep `threshold_2` as the Schmitt-style reference family and decide which current outputs or guards are genuinely Schmitt-specific rather than threshold-generic.
4. If one ambiguity remains after those decisions, use a very small trial artifact instead of widening production code immediately.

## Code-Facing Checklist

- `clode/observers/_definitions.py`, `clode/observers/types.py`, `clode/observers/metadata.py`: preserve the separate trigger-variable and feature-variable roles while shaping the next event-family declaration layer
- `clode/_opencl/source_builder.py`, `clode/_opencl/models.py`, `clode/_opencl/executors.py`: document how build signature, source preamble, and runtime-setting uploads currently divide responsibility
- `clode/kernels/observers/observer_threshold_1.clh`, `clode/kernels/observers/observer_local_maximum.clh`, `clode/kernels/observers/observer_threshold_2.clh`: use the unwired threshold stub and the live deterministic families as the immediate evidence base, but treat `observer_threshold_1.clh` as a historical sketch rather than the target layout
- `.design/reference/observer_concept_audit.md`, `.design/ideas.md`, `.design/development_roadmap.md`, `.design/package_state.md`: keep the planning and current-state narrative aligned with the threshold-topology and threshold-parameterization split

## Acceptance Criteria

- the design docs and public docs distinguish threshold crossing from Schmitt triggering as separate user-facing concepts
- the docs make clear that `threshold_2` currently uses warmup-derived fractional thresholds and Schmitt-style hysteresis, with slope gates as a narrower noise-oriented refinement
- the next implementation candidate is named with evidence: a lean two-pass fractional directional-threshold family rather than hysteresis added directly to `threshold_1`
- the design docs make clear which threshold-family follow-ons still fit the lean fixed schema
- the docs record whether mirrored `localmin` belongs in a future extremum family and whether `nhood1` is a keep-or-retire case
- generalized scalar trigger functions or hyperplane crossings remain explicitly deferred rather than drifting back into near-term scope

## Follow-on If This Lands Cleanly

1. implement the lean warmup-derived fractional directional-threshold family
2. decide when threshold-conditioned readouts become rich enough to require declaration-level schema control, and which of those belong only to Schmitt-style families
3. revisit generalized hyperplane triggers only after the landed scalar-crossing family and its next adjacent extension show what reusable trigger packaging is still missing
