# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer declaration follow-on audit

## Assumed Repo State

- Canonical public homes are `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`, with `clode._opencl` as the internal execution layer.
- The summary-observer proof slice is now landed: `Observer.summary` plus `SummaryObserverSelection` resolve to build-specialized summary variants, while `basic` and `basicall` remain compatibility presets over that same family.
- The simulation state and output ownership pass is wrapped tightly enough that integration settings, trajectory output policy, observer runtime settings, event-output policy, solver state, persistent observer state, and fetched outputs have clearer named homes on the Python side.
- Public `SolverParams`, `ObserverParams`, and the public `Stepper` enum remain thin compatibility surfaces; broader public config redesign is still deferred.
- The heavier event observers still use monolithic built-ins with large all-variable layouts, fixed feature bundles, and event-retention policy tied tightly to the chosen observer mode.

## Why this should be next

The summary slice answered the first implementation question, but it exposed the next architectural one: which observer choices belong in the build-specialized declaration path and which should remain scalar runtime settings.

That boundary matters more than immediately adding another observer because the heavier event families are not structured like summary reducers. If clODE locks in one ad hoc specialization pattern too early, it risks either macro-heavy duplication or a premature code-generation layer without knowing which parts of the observer model truly vary together.

The next useful move is therefore to audit the landed summary path, compare it to one-pass and two-pass event observers, and decide what the next actual implementation slice should be with that evidence in hand.

## Scope

- audit the landed summary-observer runtime and build path and record what part of it is a reusable observer-family pattern
- make the build-specialized versus runtime-configurable boundary explicit in the design docs
- compare the landed summary family to one-pass and two-pass event observers and identify where the same pattern does or does not generalize cleanly
- decide whether the current `#define` plus injected-preamble model is still the right near-term default, whether limited code generation is warranted anywhere, or whether some future selections should stay runtime-only
- keep actual event-observer implementation, trajectory variable-subset work, and any public custom-observer DSL out of this PR

## Key Questions

- Which observer choices necessarily change persistent-state layout or feature schema and therefore belong in the resolved spec and build key?
- Which choices are scalar enough to stay in `ObserverRuntimeSettings` or a future narrower runtime-settings object?
- Can heavier observers reuse the same declaration pattern while keeping trigger semantics, running summaries, and retained event-output policy separate enough to avoid one giant monolithic selector?
- Does the current preprocessor-based source specialization remain the best fit for the next family, or has the summary slice revealed enough repeated boilerplate to justify limited code generation?
- What is the smallest next implementation slice once those answers are written down?

## Design Constraints

- no broad public config redesign in this PR
- preserve the landed compensated time-base, accepted-step-width plumbing, solver-owned status boundary, continuation semantics, and current proof layer
- do not reintroduce solver-owned step, status, or time diagnostics through observer outputs or observer-private state
- keep semantic ownership in Python first and treat `_opencl` as the execution consumer rather than the authoritative definition layer
- no user-facing custom-observer DSL, broad inheritance hierarchy, or full code-conversion surface in the same PR
- avoid turning this PR into a benchmark, register-pressure, or "add every missing observer" campaign; the point is to record the right observer-family design boundary before more implementation
- preserve the landed summary-family behavior and current public presets while evaluating follow-on directions
- keep public docs and design notes explicit about what is clarified here and what remains deferred

## Non-goals

- no implicit or IMEX solver work
- no multi-device work
- no broader diverged-time continuation-policy redesign or matched device-side current-time model in the same PR
- no public bundled solver-stats object in the same PR
- no solver-family or stepper-model redesign in the same PR
- no public observer-state fetch API or public rename of `ObserverOutput` in the same PR
- no trajectory variable-subset implementation in the same PR
- no production event-observer implementation beyond whatever doc or prototype evidence is needed to choose the next slice
- no citation metadata, release-tag, or broader repo-surface cleanup in the same PR

## Suggested Work Slices

1. Update `.design/reference/observer_concept_audit.md` with what the landed summary slice taught about build-specialized declarations, runtime settings, and preprocessor specialization.
2. Fill and maintain a current built-in observer matrix by trigger geometry, pass structure, parameter source, retained sparse outputs, and summary-bundle scope.
3. Compare that pattern explicitly to `localmax`, neighborhood, and `threshold_2` so the next implementation slice is chosen against current code rather than intuition.
4. If one design ambiguity remains after the audit, use a very small trial artifact or scratch prototype to answer that question instead of broadening production code.

## Code-Facing Checklist

- `clode/observers/_definitions.py`, `clode/observers/types.py`, `clode/observers/metadata.py`: treat the landed summary-family declaration as the current evidence base for what belongs in a family-specific declaration layer
- `clode/_opencl/source_builder.py`, `clode/_opencl/models.py`, `clode/_opencl/executors.py`: document how build signature, source preamble, and runtime-setting uploads currently divide responsibility
- `clode/kernels/observers/observer_summary.clh`, `clode/kernels/observers/observer_local_maximum.clh`, `clode/kernels/observers/observer_threshold_2.clh`: compare the simple summary pattern to one-pass and two-pass event observers before proposing more implementation
- `.design/reference/observer_concept_audit.md`, `.design/ideas.md`, `.design/development_roadmap.md`, `.design/package_state.md`: keep the planning and current-state narrative aligned with the landed summary slice and the new audit target

## Acceptance Criteria

- the design docs explain, using current code, why the landed summary selection path is build-specialized and not just a runtime knob
- the reusable observer-family pattern from the summary slice is recorded explicitly
- the audit includes a current built-in observer matrix that distinguishes trigger geometry, pass structure, and parameter source
- the audit makes clear which parts of that pattern do not generalize directly to one-pass and two-pass event observers
- the next actual observer implementation candidate is either named with evidence or explicitly deferred behind one remaining design question

## Follow-on If This Lands Cleanly

1. pick one event-observer family whose trigger, running-summary, and retained-output policy can be split cleanly enough to test the declaration model beyond summary reducers
2. evaluate count-only versus timestamp-retaining event policies as the next likely state-footprint win for heavier observers
3. revisit limited code generation only if the next family repeats the same build-specialization boilerplate rather than because the summary slice alone felt verbose
