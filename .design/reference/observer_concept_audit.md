# Observer Concept Audit

Purpose: audit what clODE observers currently mean in code and docs, compare that model to standard ODE-package event functions, and record the current observer-state and observer-output design questions.
Read when: planning observer refactors, comparing clODE to `solve_ivp`-style events, or deciding how Python should model event detection, running feature reduction, and observer outputs.
Update when: the observer concept changes materially, a new observer-family declaration slice lands, or trajectory/output policy grows close enough to observers that the boundary needs to be restated.

## Bottom line

- clODE observers are not just SciPy-style event functions. They are compile-time-selected, stateful, on-device feature pipelines that may detect events, accumulate online statistics, retain sparse event outputs, and emit a structured readout without storing full trajectories.
- The semantic observer catalog is now much clearer than it was before the threshold-family cleanup: summary-only reducers, lean threshold and Schmitt families, `local_extremum`, and `neighborhood_return` all have explicit semantic homes, while `threshold_2`, `local_max`, and `neighbourhood_2` remain retained heavier legacy workflows.
- The current declaration boundary is also clearer: build-specialized observer choices belong in the resolved observer spec, while scalar runtime knobs stay in runtime settings when they do not change layout or feature schema.
- The heavier event observers still bundle too many concerns under one built-in mode: trigger semantics, warmup or two-pass behavior, persistent state layout, feature-output schema, and retained event-output policy.
- The active design question is no longer which threshold-family corner to add next. It is whether a shared accepted-step solution-buffer concept should precede observer bundle work, or whether the current family-local buffers are already the right level of duplication.
- `neighbourhood_1` still looks more like a keep-or-retire case than a good template for future observer generalization unless a more deterministic workflow emerges.

## Fast path

Stop after this section unless the task clearly needs the deeper observer matrix or older design argument.

- For the active planning question, read `### Accepted-step solution-buffer concept`, `### Current kernel-utilities audit`, `## Recommended next planning pass`, and `### Remaining trigger-quality gaps after the current pass`.
- For kernel observer cleanup, read `### Preferred method roles inside one kernel observer` plus `### Sampled detection versus refined event outputs`.
- For threshold, Schmitt, or neighborhood-return semantics, jump to `## Current trigger and bundle design axes`.
- For compatibility-only questions, use `.design/reference/compatibility_boundary_audit.md` instead of this longer note.
- For current delivery scope, stop and use `.design/next_pr.md` plus `.design/ideas.md`; this file is the deeper rationale layer.

## Current live model

### Python surface

- `clode.observers.types.Observer` still selects one built-in family or compatibility preset.
- The threshold families now expose semantic public names at that surface: `Observer.threshold_crossing` for the current absolute one-boundary threshold observer, `Observer.normalized_threshold_crossing` for the lean warmup-derived one-boundary threshold observer, `Observer.schmitt_trigger` for the absolute hysteresis family, and `Observer.normalized_schmitt_trigger` for the lean warmup-derived hysteresis family. `Observer.threshold_2` remains a retained legacy public name for the heavier fully featured normalized Schmitt path, while `thresh1` and `thresh3` are internal implementation names only.
- `Observer.local_extremum` and `Observer.neighborhood_return` are now also semantic public observer names, while `Observer.local_max` and `Observer.neighbourhood_2` remain retained heavier legacy workflows around the same general trigger families.
- `ThresholdCrossingConfig`, `SchmittTriggerConfig`, `LocalExtremumConfig`, and `NeighborhoodReturnConfig` are now the observer-family-specific semantic config objects beyond the summary selection surface. The threshold and Schmitt config classes are shared across observer pairs, while the extremum and neighborhood-return configs each resolve to one semantic observer.
- `Observer.summary` plus `SummaryObserverSelection` are now the first explicit family-specific declaration surface for observer feature selection.
- `ObserverParams` remains the public compatibility bundle, while `ObserverRuntimeSettings` and `EventOutputSettings` carry the narrower internal runtime settings and retained-event policy.
- The split between `e_var_ix` and `f_var_ix` is semantically important, not historical residue: users often need to trigger on one variable while measuring features on another, such as a slow calcium-like variable for event geometry and a faster voltage-like variable for amplitudes or peak counts.
- `ObserverDefinition` in `clode/observers/_definitions.py` still declares a build define, a two-pass flag, feature-name generation, and a matched persistent or event layout.
- The summary slice added a first family-specific resolved declaration path, but the definition layer still does not model event condition, running reducer family, variable-selection intent, or retained event-output intent as reusable concepts across the heavier event families.
- `FeatureSimulator` now accepts family-specific `observer_configuration` values for the threshold, Schmitt, extremum, and neighborhood-return families, while the broad `observer_*` compatibility arguments remain live.

### Runtime and build surface

- `SourceBuilder` still emits exactly one observer define for a features program build.
- `OpenCLFeatureExecutor` still owns one resolved observer spec, one observer-runtime-settings buffer, one observer-state buffer, and one feature buffer for that build.
- The landed summary slice extended `ResolvedObserverSpec` with a build variant and injected source preamble, and `BuildKey` now distinguishes observer variants beyond just the observer define plus `N_STORE_EVENTS`.
- `ObserverOutput` remains a string-keyed structured-array readout whose schema is fixed by the resolved observer spec.

### Kernel surface

- `clode/kernels/observers.cl` still expects one `ObserverState` type and one fixed function family: initialize, optional warmup, update, event test, event-feature computation, finalization, and continuation cleanup.
- The summary slice moved `basic` and `basicall` behind one `observer_summary.clh` template whose array extents and output schema are specialized by injected preamble macros.
- `local_extremum`, `localmax`, and `nhood1` remain one-pass event detectors.
- `neighborhood_return`, `thresh2`, and `nhood2` remain two-pass event detectors with a warmup solve.
- Heavy event observers still allocate state arrays scaled by `N_VAR`, `N_AUX`, and sometimes `N_STORE_EVENTS`.

### Preferred method roles inside one kernel observer

- `initializeObserverState(...)` should seed persistent buffers, counters, and retained event-output storage.
- `warmupObserverState(...)` should collect only the warmup-pass geometry needed to parameterize the live detector.
- `initializeEventDetector(...)` should convert that warmup geometry plus the rewound initial live-pass sample into concrete detector state.
- `updateObserverState(...)` should advance accepted-step history buffers and continuous per-step reducers before the event test runs.
- `eventFunction(...)` should be the primary place that decides whether the observer's public event fired on the current step. If a refined timestamp or exit fraction is needed, that trigger geometry should stay local to this phase or to `computeEventFeatures(...)` rather than turning `updateObserverState(...)` into a proxy event detector.
- `computeEventFeatures(...)` should consume the already-available event-local geometry and update sparse outputs, counts, periods, or terminal-event state.
- `finalizeFeatures(...)` and `finalizeObserverState(...)` remain the readout and continuation cleanup phases.

That preferred split now matches the lean threshold families, `local_extremum`, the cleaned-up `neighborhood_return`, and the updated `nhood2` exit path more closely than the earlier pending-event workaround did.

### Sampled detection versus refined event outputs

- clODE should distinguish event detection from event-output refinement.
- Detection decides whether the observer's public event fired on the accepted-step sampled history, for example from a sign change, threshold bracket, or state-machine transition.
- Refinement happens only after that sampled event is accepted. It can improve stored metrics such as event time, event value, or elapsed-time period endpoints by interpolating over the same bracketing samples.
- Refinement should not silently redefine the detection contract. A family that wants true root-finding-based detection rather than sampled bracketing should state that explicitly as a different semantic model.
- This distinction explains the current intended model for threshold and neighborhood-return observers: sampled detection on the live step history, then interpolated timestamps for stored event outputs.

### Accepted-step solution-buffer concept

- The minimal shared solution-buffer concept is a `K`-step history of accepted samples `(t, x[N_VAR], dx[N_VAR])` containing the current sample and the previous `K - 1` samples.
- Some observers also need mirrored solve-relative elapsed-time buffers or auxiliary-value buffers, but those are extensions on top of the core accepted-step history rather than the core abstraction itself.
- `K` is family-specific. Threshold and neighborhood-exit refinement usually need two samples, while extremum refinement usually needs three.
- That accepted-step buffer is a compact view of the recent local trajectory geometry. It is enough for many event-detection and storage-refinement tasks, but it is not yet the same thing as general dense output.
- A future shared solution-buffer abstraction should describe this accepted-step geometry once and let observer families consume the slice size they actually need.

### Helper-utility placement guidance

- Prefer `clode/kernels/clODE_utilities.cl` for reusable math or geometry primitives when the operation is family-agnostic and its address-space contract can stay explicit and portable.
- Passing explicit sizes as helper arguments, as in `norm_inf(..., N_VAR)`, is a good default when it keeps the helper reusable across steppers and observers.
- Keep helper code local to an observer CLH when it bakes in observer-specific policy, output semantics, or awkward coupling to one concrete `ObserverState` layout.
- OpenCL address-space rules are part of the helper API surface. A candidate shared helper that cannot be expressed cleanly without layout-specific or address-space-specific assumptions may be better left local until a generic signature is clear.
- The design target is a small shared utility layer for common array math, interpolation primitives, and eventually solution-buffer updates, plus thinner observer-local policy code rather than a utilities file full of observer-specific semantics.

### OpenCL address-space guidance for helpers

- Work-item-local arrays and struct members in clODE observers, steppers, and RNG state are `__private` objects.
- The OpenCL spec makes `__generic` an optional feature: it is baseline only for OpenCL C 2.0 and otherwise requires the OpenCL C 3.0 generic-address-space feature path. clODE should not rely on `__generic` as the default portability contract for shared helpers.
- Local build probes on the current NVIDIA OpenCL path were enough to establish the safe house style: explicit `__private` helpers work with both stack arrays and struct-member arrays, and explicit `__global` helpers work for buffer-backed arrays.
- Unqualified pointer or array parameters may compile on some toolchains, but they are a weak portability signal and should not be the intended contract for new shared helpers.
- House style going forward: scalar-only helpers are always the easiest to share; pointer or array helpers should name the expected address space when they cross file boundaries; if the same algorithm is needed in multiple named address spaces, prefer thin wrappers or caller-side private scratch copies over hidden dependence on generic-pointer behavior.

### Current kernel-utilities audit

- Shared scalar helpers are already in the right place: compensated time and integral helpers, running means, and scalar interpolation routines are address-space-neutral and broadly reusable across steppers and observers.
- Current array helpers in `clODE_utilities.cl` such as `norm_*`, `array_*`, `quadraticInterpVertex*`, and the three-sample extremum helpers are effectively caller-private utilities today. Their live call sites feed them stack arrays, step-local scratch buffers, or observer-local temporary copies rather than raw non-private buffers.
- That current usage is a useful rule, not just an accident: if a helper needs contiguous array access but the source data lives inside a larger observer-state layout, a caller-side private slice is often the cleanest contract until a broader solution-buffer abstraction lands.
- Observer-family geometry that still depends on one concrete `ObserverState` layout or one family's event semantics, such as the current neighborhood-return normalized-ball exit calculations, should remain local to the observer CLH for now.
- The next utilities pass should therefore focus on two low-risk buckets: explicit address-space policy for shared array helpers, and shared accepted-step solution-buffer update primitives once the observer-side `K`-sample buffer shape is promoted into a common abstraction.

## Current built-in observer matrix

Pass count alone is not enough to classify the current observers. The code currently exposes three parameter-source classes:

- no learned internal parameters at all
- on-the-fly internal-parameter discovery during the live pass
- warmup-derived internal parameters from a first pass

That distinction matters because it changes what can stay in runtime settings, what likely belongs in a family-specific declaration, and what makes a trigger geometry stable or unstable over the solve window.

| Public modes | Trigger class | Pass structure | Internal-parameter source | First warmup pass required? | Retained sparse output | Persistent summary scope today | Design notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Observer.summary`, `Observer.basic`, `Observer.basic_all_variables` | none; online reduction only | one pass | none beyond the resolved declaration and scalar runtime selectors | no | none | selected summary subset for `summary`; preset subset for `basic`; all-state plus all-aux preset for `basic_all_variables` | cleanest family; selection changes layout and schema directly |
| `Observer.local_extremum` | local extremum in the selected variable via derivative sign change with three-sample refinement | one pass | static runtime selectors only; polarity stays in runtime settings | no | local-extremum timestamps and values up to `N_STORE_EVENTS` | lean event stream only | landed semantic extremum family; confirms polarity can stay runtime-side when the readout schema stays fixed |
| `Observer.local_max` | local maximum in `fVarIx` via derivative sign change with three-sample refinement | one pass | static runtime selectors only | no | local max/min timestamps and values up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle | retained heavier maxima-oriented workflow when IMI, amplitude, and summary outputs are still needed |
| `Observer.neighbourhood_1` | entry into an `nHoodRadius` ball around `x0` | one pass | `x0` discovered online from the first local minimum in `eVarIx`; normalization range keeps evolving during the same pass | no formal warmup, but yes to live-pass parameter discovery | none beyond event count | all-state plus all-aux summary bundle | awkward family: trigger geometry depends on mutable normalization data during the live pass; weak evidence for keeping it as a future template |
| `Observer.neighborhood_return` | exit from an `nHoodRadius` ball around `x0` in normalized state space | two pass | warmup-derived trajectory ranges and threshold; `x0` is then pinned in the live pass as the first threshold-qualified point on `eVarIx` | yes | exit timestamps up to `N_STORE_EVENTS` | lean event stream only | landed semantic neighborhood-return family; keeps the specialized periodicity trigger without the broader legacy readout bundle, now with sampled anchors plus interpolated exit times |
| `Observer.neighbourhood_2` | exit from an `nHoodRadius` ball around `x0` in normalized state space | two pass | warmup-derived trajectory ranges and threshold; `x0` is then pinned in the live pass as the first threshold-qualified point on `eVarIx` | yes | exit timestamps up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle, plus range and `x0` readout | retained heavier periodicity-oriented workflow when the broader period, peak, and summary outputs are still needed; now shares sampled-anchor plus interpolated exit timing with `Observer.neighborhood_return` |
| `Observer.threshold_crossing` | absolute threshold crossing on `eVarIx` with runtime direction selection | one pass | static runtime selectors only | no | one threshold-event timestamp stream up to `N_STORE_EVENTS` | lean event stream only; no all-state summary bundle | confirms that a fixed-schema event family can keep threshold value and direction in runtime settings |
| `Observer.normalized_threshold_crossing` | warmup-derived normalized threshold crossing on `eVarIx` with runtime direction selection | two pass | warmup-derived fraction of the global `x` range on `eVarIx` | yes | one threshold-event timestamp stream up to `N_STORE_EVENTS` | lean event stream only; no all-state summary bundle | useful when oscillation amplitudes drift across a sweep but the workflow still wants one directional threshold stream rather than Schmitt-state outputs |
| `Observer.schmitt_trigger` | absolute up/down thresholds on `eVarIx` with Schmitt-style up/down state | one pass | static runtime selectors only | no | up/down transition timestamps up to `N_STORE_EVENTS` | lean up/down event stream plus event count | useful when the hysteresis band has a stable physical scale and only trigger timing is needed |
| `Observer.normalized_schmitt_trigger` | warmup-derived normalized up/down thresholds on `eVarIx` with Schmitt-style up/down state | two pass | warmup-derived fractions of global `x` and `dx` ranges on `eVarIx` | yes | up/down transition timestamps up to `N_STORE_EVENTS` | lean up/down event stream plus event count | useful when oscillation amplitudes drift across a sweep but the workflow still wants a lean Schmitt-style readout |
| `Observer.threshold_2` | warmup-derived normalized up/down thresholds on `eVarIx` with Schmitt-style up/down state; local maxima in `fVarIx` contribute derived summaries | two pass | warmup-derived fractions of global `x` and `dx` ranges on `eVarIx` | yes | up/down transition timestamps up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle, plus period, duty, and active-dip bundle | retained legacy full workflow when the broader oscillation-oriented readout bundle is still needed |

## Current trigger and bundle design axes

### Axes that still matter

- Trigger topology remains a first-class user-facing distinction: threshold crossing, Schmitt hysteresis, local extrema, and neighborhood return answer different questions and should not be collapsed into one conceptual ladder.
- Parameter source matters alongside topology: static runtime thresholds, warmup-derived fractions, and live-pass discovery create different stability and packaging needs even when the user-facing trigger looks similar.
- Readout weight is separate from trigger semantics: the lean semantic families and the retained heavy legacy families now coexist, and the next bundle question is about shared accepted-step buffers and recurring readouts rather than about relitigating the threshold-family split.
- Variable roles stay explicit. `eVarIx` and `fVarIx` should remain distinct because real workflows often trigger on one variable while measuring period, amplitude, or extrema on another.

### Current family guidance

- Keep threshold crossing and Schmitt triggering as separate user-facing concepts even when internal helpers overlap.
- Keep family-specific semantic config objects as the preferred surface, with `ObserverParams` and the `observer_*` keywords remaining compatibility adapters.
- Treat `threshold_2`, `local_max`, and `neighbourhood_2` as retained legacy evidence for bundle decisions, not as the default template for new semantic families.
- Keep `neighbourhood_1` under explicit keep-or-retire review rather than using it to drive future abstraction.
- The archived threshold-family rollout argument and earlier family-ranking analysis now live in `../archived/reference_surface_refresh_2026_05_28/observer_threshold_family_background.md` so this live note can stay focused on current semantics and follow-on questions.

## Relation to `solve_ivp` events

`scipy.integrate.solve_ivp` events are scalar root functions of `(t, y)` with optional `direction` and `terminal` semantics. The solver finds zero crossings, returns `t_events` and `y_events`, and otherwise leaves summary-feature logic to user code layered on top of the solver.

That model overlaps with clODE only at the event-trigger layer:

- `threshold_2`, `local_max`, and the neighborhood observers have an event-condition component that is broadly comparable to `solve_ivp` events.
- `basic`, `basicall`, and `summary` do not correspond to event functions at all; they are online reducers over the whole solve window.
- clODE observers additionally own persistent per-trajectory state, online summary updates, event-count or timestamp retention, and a final structured readout. Those are central to clODE's large-ensemble workflow and should not be reduced to plain callback-style events.

Useful semantics to borrow from `solve_ivp` as naming or modeling guidance:

- explicit distinction between event condition and event readout
- explicit directionality for crossings
- explicit terminal versus non-terminal event behavior
- clear statement that multiple crossings inside one coarse step may still be missed unless the kernel tracks more geometry
- a scalar trigger function or hyperplane crossing as the natural long-term generalization of threshold-style events

Semantics that should not become clODE's primary observer model:

- arbitrary Python callbacks inside the timestep loop
- dense output as the defining event abstraction
- treating all feature extraction as just event detection

Near-term design consequence:

- generalized scalar trigger functions are still a good long-term umbrella, but the narrow scalar-crossing proof is already landed; any broader trigger-function DSL or hyperplane generalization should remain later work after the current solution-buffer and bundle questions are settled.

## What the landed summary slice taught

### Why custom summary selections are build-specialized

- `SummaryObserverSelection` does not just rename outputs. It changes the ordered feature schema, the persistent `ObserverState` layout, and the generated struct name used for PyOpenCL dtype matching.
- Because those choices change both kernel-visible array extents and host-visible dtype layout, they currently have to participate in the resolved observer spec and the program-cache key.
- In practice, “build-specialized” here means the Python declaration resolves to one concrete summary variant, `SourceBuilder` injects a preamble that defines the required counts and index arrays, and `OpenCLFeatureExecutor` rebuilds when that resolved variant changes.
- The runtime-settings buffer remains the right home for scalar thresholds, counters, and distinguished variable indices that do not alter feature schema or persistent-state shape.

### Systematic pattern worth reusing

The summary slice exposed a reusable five-part pattern:

1. one Python-side declaration object for the family-specific selection surface
2. one resolved spec that derives feature names, persistent layout, and build identity together
3. one build-signature path for anything that changes kernel-visible shape or output schema
4. one narrower runtime-settings path for scalar controls that do not change layout
5. one executor invalidation rule that distinguishes rebuilds from simple runtime-setting uploads

That pattern is likely reusable across observer families even if the concrete declaration objects differ.

### Where that pattern stops being simple

- Summary reducers are unusually clean because the update loop is separable: selected reductions can be expressed as independent per-step accumulators with no warmup pass and no event-retention geometry.
- `localmax`, neighborhood, and threshold observers couple trigger semantics, buffered geometry, running period or amplitude summaries, and retained sparse event outputs much more tightly.
- For those heavier families, “selection” likely needs to be split into separate concepts such as event trigger, tracked running-summary groups, and retained event-output policy rather than copied directly from `SummaryObserverSelection`.
- The summary-family preamble is also verbose because it uses `#define` counts and compile-time index arrays. That is acceptable for one family, but it would become hard to maintain if every event family duplicated the same pattern ad hoc.

## Remaining friction points

### 1. Event-observer identity is still monolithic

- One event-observer name still chooses trigger semantics, persistent state layout, readout schema, warmup policy, and build-time kernel selection all at once.
- That is manageable for a small built-in catalog, but it blocks finer-grained feature selection and makes built-in extension feel heavier than it needs to be.

### 2. Runtime settings are broad and partially inactive

- `ObserverRuntimeSettings` still exposes one shared set of indices, thresholds, and counters for all built-ins.
- Many fields are irrelevant for a given observer, so the Python surface still does not tell users which knobs are structural and which are inert for the current mode.

### 3. Feature-selection granularity is still coarse for event observers

- The landed summary family now supports selected state, auxiliary, and slope subsets through `SummaryObserverSelection`.
- The heavier event observers still carry all-state summary arrays, even when a user may only care about period, count, or extrema on one or two variables.
- There is still no current way to ask for “count local extrema without storing the full all-variable summary bundle” or similar narrower event-observer readouts through one reusable Python definition model.

### 4. Memory footprint still follows observer mode more than requested readout

- The landed summary family now shrinks persistent state with the selected summary subset.
- `localmax`, `nhood1`, `nhood2`, and `thresh2` still size persistent arrays directly from `N_VAR` and often `N_AUX`.
- Event-capable observers also scale state or retained outputs with `N_STORE_EVENTS`.
- Today the main user control over event-observer footprint is still indirect: pick a lighter observer, reduce `max_event_timestamps`, or avoid auxiliary-heavy models.

### 5. Output schema is still fixed and stringly typed

- `ObserverOutput` is useful for inspection, but the schema is still just an ordered list of strings paired with one flat structured array.
- That makes it easy to expose results, but harder to declare optional feature groups, retained-event policies, or memory-affecting output choices on the Python side before the build happens.

## `#define` flags versus generated source versus runtime configuration

- The current `#define` plus injected-preamble model is defensible when a family has a small, testable kernel template and the selection directly changes array extents or output schema.
- A pure runtime-configuration model is only a good fit when the shape is fixed and the choice only changes scalar thresholds, selected event variables, or behavior inside an already-allocated template.
- Full source code generation is not yet justified by the current evidence. The main problem today is not that the summary slice used preprocessor specialization, but that the heavier families still lack a comparable declaration model to tell us which choices actually need specialization.
- The most plausible near-term stance is hybrid: keep Python-owned declaration objects, resolved specs, and one compiled observer family per build; keep runtime buffers for scalar knobs; consider limited code generation only if multiple observer families end up repeating the same preamble boilerplate after the next audit.

## Recommended next planning pass

- The deterministic static-trigger slice and its adjacent follow-ons are now landed: absolute and warmup-derived threshold crossing plus absolute and warmup-derived Schmitt triggering are all live behind shared semantic config classes.
- The immediate threshold-design question is no longer how to expose the semantic split or which threshold corner to add next.
- The compatibility boundary is now documented. The remaining planning question is not another naming cleanup.
- The semantic Schmitt families and the retained `threshold_2` timing path now route both up and down transitions through `eventFunction(...)` plus `computeEventFeatures(...)`, so the clearest resolved method-role asymmetry is no longer in that slice.
- Utilities currently look mostly like private-use implementation support, not the right standalone next PR.
- The next concrete planning target should be deciding whether a shared accepted-step `K`-sample solution-buffer concept is the missing prerequisite for clearer observer-state and observer-output bundle work.
- That planning pass should say explicitly whether the oscillation-oriented readout seam should follow immediately after that audit, or whether family-local buffers are already clear enough that bundle work can proceed without a shared buffer abstraction.
- `localmax` should be reconsidered as one half of a likely extremum family with a max/min polarity switch.
- `nhood2` should remain documented as a specialized periodicity detector for complex limit cycles even though it is not the next proof target.
- `nhood1` should remain under explicit keep-or-retire review and should not define the future observer abstraction unless a clearer deterministic workflow emerges.
- Generalized `solve_ivp`-style scalar trigger functions or hyperplane crossings should remain a later generalization after the landed scalar-crossing slice and the next adjacent follow-on clarify what reusable trigger packaging is still missing.
- Trajectory variable-subset storage should remain separate output-policy work and stays tracked in `.design/ideas.md`.

## Event-function audit follow-up

### Findings from the current semantic kernels

- The semantic Schmitt kernels now follow an explicit x-only contract. `Observer.schmitt_trigger` and `Observer.normalized_schmitt_trigger` both treat `x_down_threshold` as a literal configured boundary, so zero is representable, equality with `x_up_threshold` is allowed, and there is no derivative-gate sentinel behavior on the semantic path.
- Derivative gates now live only on the retained legacy `Observer.threshold_2` compatibility surface. That keeps the semantic Schmitt family predictable while preserving the older noisy-trace workflow where it already existed.
- The semantic Schmitt families and the retained `Observer.threshold_2` path now follow the preferred event-method split more closely: `updateObserverState(...)` advances accepted-step history and continuous reducers, `eventFunction(...)` detects both state-machine transitions on sampled history, and `computeEventFeatures(...)` refines and stores the transition-specific outputs.
- `thresholdTransitionTime(...)` now returns the later active-gate time even when one threshold gate was already active at the start of the bracketing step. That keeps the retained slope-gated `threshold_2` path consistent with the intended "wait for the second active gate" semantics.
- The semantic config layer now validates the main trigger geometry that used to be implicit. Schmitt rejects `x_up_threshold < x_down_threshold`, normalized threshold-like families reject fractions outside `[0, 1]`, and `NeighborhoodReturnConfig.radius` must be positive.
- `min_amp` is currently not one unified semantic concept across the lean trigger families. In the one-pass absolute threshold and absolute Schmitt observers it acts as a live-pass range gate on the event variable, while in the normalized families the same field is compared against a warmup-derived full-window amplitude. That difference is defensible, but it is not obvious from the current config names alone.
- `Observer.local_extremum` is currently the clearest and most predictable lean event family. Its trigger is a sampled derivative sign change on the selected variable, and its timestamps and values are then refined with the bounded three-sample quadratic helper. The main caveat is that it still detects sampled sign changes, not arbitrary zero roots of the derivative.
- `Observer.neighborhood_return` is still the least intuitive semantic family even though the core trigger is better behaved than legacy `nhood1`. The anchor point `x0` is still latched to the first sampled point below the anchor threshold rather than to an interpolated threshold hit, zero-range dimensions are still skipped in the normalized distance, and the new exit interpolation now follows the linearly interpolated full-state segment between the last inside sample and the first outside sample of the radius ball.
- `Observer.neighbourhood_2` now shares that improved exit interpolation, including interpolated elapsed-time periods, but it still inherits the heavier legacy state bundle and some older step-buffer conventions.
- The current tests now cover equal-threshold Schmitt configuration, inverted Schmitt rejection, semantic-Schmitt dx rejection on the compatibility path, direct lean-Schmitt up/down transition storage, later-active-gate timing on the retained `threshold_2` path, out-of-range normalized thresholds, positive-radius neighborhood validation, and interpolated exit timing for `Observer.neighborhood_return` and `Observer.neighbourhood_2`. Exact-boundary conventions, zero-range warmup behavior, and the sampled-anchor semantics of the neighborhood-return family still need more coverage.

### Remaining trigger-quality gaps after the current pass

- The clearest remaining structural observer gap is no longer the Schmitt method-role split. It is that several observers still carry family-local accepted-step buffers and near-duplicate buffer-update patterns, which makes later observer-state and observer-output bundle work harder to compare than it needs to be.
- `Observer.local_extremum` and `Observer.local_max` now have good timestamp refinement, but they still detect sampled derivative sign changes rather than solving for derivative roots directly.
- `Observer.neighborhood_return` and `Observer.neighbourhood_2` now have better exit timestamps, but their anchors remain sampled rather than interpolated, and the exact equality conventions around the anchor threshold and neighborhood radius still deserve tighter tests.
- `Observer.neighbourhood_1` still looks like the weakest live trigger family because its normalization geometry evolves during the same pass that uses it.

### Recommended improvements before public figure-heavy docs

1. Add edge-case tests for the remaining semantic ambiguities before widening the docs story further: exact-threshold boundary cases, zero-range warmup dimensions, and the precise anchor-and-exit timing of the neighborhood-return family.
2. Keep the new visualization examples aligned with the code path they explain: threshold interpolation for the threshold families, x-only state-machine entry and exit for semantic Schmitt, three-sample extremum refinement for `local_extremum`, and sampled-anchor plus normalized-ball exit interpolation for `neighborhood_return` and `nhood2`.
3. Keep the public docs precise about which event times are interpolated and which are sampled. Threshold and semantic Schmitt times are refined within the step, `threshold_2` still has the broader legacy slope-gated path, `local_extremum` uses three-sample refinement, and the neighborhood-return families now keep sampled anchors but refine exit times within the step.

## Audit anchors

- `clode/observers/_definitions.py`
- `clode/observers/types.py`
- `clode/observers/metadata.py`
- `clode/simulation/features.py`
- `clode/simulation/results.py`
- `clode/_opencl/observer_metadata.py`
- `clode/_opencl/executors.py`
- `clode/_opencl/source_builder.py`
- `clode/_opencl/models.py`
- `clode/kernels/observers.cl`
- `clode/kernels/observers/observer_summary.clh`
- `clode/kernels/observers/observer_local_maximum.clh`
- `clode/kernels/observers/observer_schmitt_trigger.clh`
- `clode/kernels/observers/observer_threshold_2.clh`
- `.design/reference/compatibility_boundary_audit.md`
- `clode/kernels/features.cl`
- `clode/kernels/initializeObserver.cl`
- `docs/feature_extraction.md`
- `scipy.integrate.solve_ivp` event documentation
