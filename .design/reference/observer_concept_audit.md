# Observer Concept Audit

Purpose: audit what clODE observers currently mean in code and docs, compare that model to standard ODE-package event functions, and record what the landed summary-observer slice taught about declaration structure, build specialization, and runtime configuration.
Read when: planning observer refactors, comparing clODE to `solve_ivp`-style events, or deciding how Python should model event detection, running feature reduction, and observer outputs.
Update when: the observer concept changes materially, a new observer-family declaration slice lands, or trajectory/output policy grows close enough to observers that the boundary needs to be restated.

## Bottom line

- clODE observers are not just SciPy-style event functions. They are compile-time-selected, stateful, on-device feature pipelines that may detect events, accumulate online statistics, retain sparse event outputs, and emit a structured readout without storing full trajectories.
- The first declaration-model slice is now landed for summary-only reducers: `Observer.summary` plus `SummaryObserverSelection` resolve to a build-specialized summary variant, while `Observer.basic` and `Observer.basic_all_variables` remain compatibility presets over that family.
- That landed slice exposed a reusable pattern: a Python-side declaration resolves to one concrete spec that determines feature names, persistent layout, and build identity together, while scalar runtime knobs remain in the runtime-settings path.
- The first deterministic static-trigger slice is now also landed: `Observer.threshold_crossing` provides one-pass absolute threshold crossings with a lean timestamp-plus-count readout and runtime direction selection.
- The threshold catalog now spans four public semantic combinations: absolute threshold crossing, warmup-derived fractional threshold crossing, absolute Schmitt hysteresis, and warmup-derived fractional Schmitt hysteresis.
- The threshold-family semantic config seam is now also landed for the current built-ins: `ThresholdCrossingConfig` is shared by the two one-boundary threshold observers, `SchmittTriggerConfig` is shared by the two Schmitt families, and `ObserverParams` plus the `observer_*` keyword arguments remain compatibility surfaces.
- `Observer.local_extremum` is now landed as the lean polarity-selectable extremum family, while `Observer.local_max` remains the retained heavier maxima-oriented workflow.
- `Observer.neighborhood_return` is now landed as the lean two-pass normalized neighborhood-return family, while `Observer.neighbourhood_2` remains the retained heavier periodicity-oriented workflow.
- The heavier event observers still bundle too many concerns under one built-in mode: trigger semantics, warmup or two-pass behavior, persistent state layout, feature-output schema, and retained event-output policy.
- That threshold slice confirmed the intended boundary from the opposite side of the summary family: scalar trigger controls such as crossing direction can stay runtime-side when the readout schema remains fixed.
- That threshold-family follow-on is now also landed: `Observer.normalized_threshold_crossing` provides a lean warmup-derived normalized directional-threshold observer with the same timestamp-plus-count readout shape as `Observer.threshold_crossing`.
- `Observer.schmitt_trigger` and `Observer.normalized_schmitt_trigger` are now landed as lean absolute and warmup-derived Schmitt families, while `Observer.threshold_2` names the retained heavier warmup-derived normalized Schmitt path.
- `neighbourhood_1` currently looks more like a keep-or-retire case than a good template for future observer generalization unless a more deterministic workflow emerges.

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

## Near-term candidate family ranking

| Candidate family | Current evidence | Near-term fit | Design read |
| --- | --- | --- | --- |
| Absolute threshold crossing in state-variable units | now landed as `Observer.threshold_crossing` with a lean event stream and runtime direction selection | landed first proof target | confirms that not every observer-family extension needs schema-changing specialization |
| Warmup-derived normalized threshold crossing | now landed as `Observer.normalized_threshold_crossing` with a lean one-stream readout and warmup-derived threshold parameterization | landed second threshold proof target | separates threshold parameterization from Schmitt topology without overloading the Schmitt family |
| Absolute Schmitt trigger | now landed as `Observer.schmitt_trigger` with a lean up/down timestamp readout and thresholds expressed directly in event-variable units | landed follow-on proof target | confirms that trigger topology and threshold parameterization can vary independently while still sharing one public config model |
| Local extremum family | now landed as `Observer.local_extremum` with `ExtremumPolarity` selection while `Observer.local_max` remains the retained heavier legacy workflow | landed adjacent proof target | confirmed that extremum polarity can stay runtime-side while the lean event schema stays fixed |
| Generalized hyperplane crossing | conceptually matches the useful part of `solve_ivp`-style scalar events | later | promising long-term trigger umbrella, but it needs a declaration for scalar trigger functions or hyperplane parameters that the current runtime path does not yet package cleanly |
| Neighborhood return triggers | now landed as `Observer.neighborhood_return` for the lean specialized path while `Observer.neighbourhood_2` remains the retained heavier legacy workflow | landed specialized follow-on | keeps `nhood2`'s better-behaved trigger semantics available without forcing the broader periodicity bundle onto the lean path; `nhood1` remains the clearest retirement candidate |

## Trigger and parameter taxonomy

### 1. Static trigger families

- Summary reducers and `local_max` do not need a first pass to discover internal trigger parameters.
- Their runtime knobs are mostly scalar selectors, thresholds, or counters.
- These families are the easiest place to separate trigger geometry from retained outputs or running summary groups.
- A one-pass absolute-threshold observer belongs in this class if its thresholds are interpreted directly in state-variable units instead of being derived from a warmup pass.

### 2. Online-discovery trigger families

- `neighbourhood_1` does not run a formal warmup pass, but it still discovers internal trigger parameters during the live pass.
- The anchor point `x0` is set from the first local minimum in `eVarIx`, and the normalization basis in `eventFunction(...)` uses trajectory ranges that are still changing over time.
- That makes `nhood1` a distinct design case rather than just another one-pass detector.

### 3. Warmup-derived trigger families

- `threshold_2` and `neighbourhood_2` both need a first pass to establish internal geometry before the main event pass starts.
- That is a good fit for a two-pass family concept, but the current code still mixes warmup-derived geometry with live-pass feature bundles and retained-output policy under one built-in name.

### Threshold topology versus threshold parameterization

- The current threshold discussion is clearer if it is split along two independent axes instead of treated as one linear family ladder.
- The first axis is trigger topology: a single-boundary directional crossing versus a two-boundary Schmitt-style state machine.
- The second axis is threshold parameterization: absolute state-variable units versus warmup-derived fractional position inside the observed event-variable range.
- `Observer.threshold_crossing` currently occupies the absolute plus single-boundary corner.
- `Observer.schmitt_trigger` currently occupies the absolute plus Schmitt-style corner.
- `Observer.normalized_schmitt_trigger` currently occupies the warmup-derived plus Schmitt-style corner.
- The warmup-derived normalized directional-threshold corner is now landed as `Observer.normalized_threshold_crossing`, with the same lean readout style as `Observer.threshold_crossing`.
- Landing the absolute Schmitt corner confirmed that trigger topology and threshold parameterization can stay separate in the public model even when the internal machinery overlaps.

### Why threshold crossing and Schmitt triggering should stay separate user-facing concepts

- Threshold crossing answers a one-boundary question: when does `x` pass this level, in which direction?
- Schmitt triggering answers a stateful two-boundary question: when does `x` enter and leave a latched upstate without counting chatter as repeated events?
- It is true that a Schmitt trigger with `xUp == xDown` and both directions collapses to a degenerate single-boundary case, but that equivalence is better treated as an implementation detail than as the public UX model.
- Separate public concepts keep the docs legible, reduce inert parameters on the simple path, and make it easier to explain why a lean threshold family and a heavier hysteresis family naturally want different default readouts.
- Internal helper sharing is still compatible with this separation. clODE can reuse trigger machinery without forcing users to think in Schmitt-state terms when they only want a directional crossing.

### Why parameter packaging now matters more than one more threshold kernel

- The current public compatibility surfaces still expose one broad observer bundle and one broad set of `observer_*` keyword arguments.
- That surface is now the most obvious UX mismatch in the threshold family: a one-boundary threshold crossing only needs one threshold scalar plus direction, while a Schmitt trigger naturally wants two value thresholds and optionally two slope thresholds.
- This is also where the current names become awkward. `x_up_threshold` is a tolerable compatibility spelling in a generic bundle, but it is not the right semantic knob name for a downward-only threshold crossing.
- The architectural implication is that clODE should move toward observer-family-specific config surfaces that show only the relevant knobs for the chosen observer while keeping `ObserverParams` and the `observer_*` constructor keywords as compatibility layers.
- That config split should be semantic rather than pass-count-driven: threshold-crossing families should expose one `threshold` plus `direction`; Schmitt-style families should expose `x_up_threshold`, `x_down_threshold`, and optional slope gates.
- That semantic config seam was enough to land the full current four-corner threshold catalog without worsening the naming debt, but it also means the selected observer rather than the config class now owns absolute versus warmup-derived interpretation.

### What the landed threshold-family config seam proved

- Family-specific config objects are enough to improve the UX without changing the executor or kernel contracts. The Python side can adapt semantic configs into the existing normalized runtime settings.
- Once one config class is shared across an absolute and a warmup-derived pair, the config type can no longer determine the observer family by itself. The selected observer remains the semantic source of truth.
- Some knobs can still recur across families, such as `min_amp`, without forcing one giant public bundle. The better rule is to surface a knob where it is active for that family, not where it might someday become reusable elsewhere.
- The current four-family threshold catalog makes the next boundary clearer: recurring controls such as `min_amp` and `max_event_count` should stay on the active oscillation-oriented family configs until output-bundle work shows a better shared seam.
- The broad `ObserverParams` bundle and `observer_*` keywords are still useful compatibility layers, but they no longer need to define the preferred semantic vocabulary for threshold families.

### Why neighborhood-return triggers remain valuable

- `nhood2` is materially better behaved than `nhood1` because the warmup pass fixes the trajectory ranges first, and the live pass only accepts `x0` after `x[eVarIx]` drops below a threshold derived from that warmup.
- That matters when `eVarIx` is a slow variable used to identify a meaningful phase of a bursting or otherwise multiscale oscillation.
- The family is therefore not just a threshold surrogate. Its intended semantics are to define a small neighborhood around a good full-state anchor point and then detect departure after a return to that neighborhood.
- That can be a more robust periodicity detector than simple threshold crossings or local-extremum triggers when complex limit cycles revisit similar scalar values at several distinct phases.
- The design conclusion is narrower than “neighborhood triggers are bad.” That read now lands as `Observer.neighborhood_return` for the lean specialized path, with `Observer.neighbourhood_2` retained as the heavier periodicity-oriented workflow and `nhood1` still the weak link that should not drive the general observer abstraction.

## Design implications from the trigger matrix

### Why parameter source is a first-class design axis

- Pass count alone does not explain what kind of declaration a family needs.
- The more useful axis is where the trigger geometry gets its internal parameters: nowhere, from static runtime knobs, from live-pass discovery, or from a warmup pass.
- That axis should likely sit alongside trigger geometry, running-summary groups, and retained sparse-output policy in the next observer design.

### Why `eVarIx` versus `fVarIx` must stay explicit

- The observer catalog should not collapse event geometry and feature readout onto one implied primary variable.
- A real workflow is to trigger on a slower or more reliable reference variable while measuring period, amplitude, maxima, or other summaries on a different faster variable.
- That means future threshold, extremum, and hyperplane-trigger families should preserve separate trigger-variable and feature-variable roles, even if they default to the same variable in simpler use cases.

### Why absolute threshold crossing was the right counterexample to summary-style specialization

- The summary family proved the shape-changing side of the declaration boundary: selected reductions change layout and schema, so they must be build-specialized.
- An absolute-threshold family is the opposite kind of test: threshold values, hysteresis values, and crossing direction are scalar behavior knobs that should usually stay in runtime settings because they do not inherently change observer-state layout.
- That makes a one-pass absolute-threshold family the best next proof target for checking whether the declaration model can stay narrow when the readout schema is fixed.
- It is also more deterministic and more intuitive to users than the current neighborhood-return ideas.

### What the landed `threshold_crossing` slice taught

- A lean event family can keep a fixed output schema while still exposing a meaningful new runtime control, here `EventDirection`.
- The first slice did not need build-signature changes beyond selecting the observer family itself: changing absolute threshold values or crossing direction does not rebuild the program.
- Rewriting the old `observer_threshold_1.clh` sketch down to one timestamp stream plus count was the right call. The monolithic all-state snapshot design was not needed to get useful behavior.
- The first landed slice still leaves one genuine follow-on question open: how much more event-conditioned readout can be added before the family stops being runtime-only and starts needing declaration-level schema control.

## What the minimum absolute-threshold family proved

The landed `threshold_crossing` slice did not clone `threshold_2`'s full readout bundle into a one-pass shell. That smaller proof target was the right call.

### Landed family shape

- Trigger geometry: scalar crossing on `eVarIx` in absolute state-variable units.
- Required runtime knobs: `eVarIx`, crossing direction, one absolute threshold, `maxEventCount`, and retained-timestamp capacity.
- Preserved runtime role: `fVarIx` remains independently selectable in the runtime contract even though the lean first slice does not yet attach a richer event-conditioned readout to it.
- Retained outputs: event count plus optional crossing timestamps, with no full copied event-state snapshot.
- Running summaries: none in the first slice; the family stays lean instead of dragging along the all-state plus all-aux bundle.

### What should stay out of the first threshold slice

- no generalized scalar trigger-function or hyperplane surface yet
- no automatic normalization from a warmup pass
- no mandatory all-state event snapshots or copied `xThisEvent`/`dxThisEvent`/`auxThisEvent` arrays
- no bundled local-extrema side channel unless a concrete readout requires it

## What `threshold_2` still teaches that `threshold_crossing` and `normalized_threshold_crossing` do not

- Warmup-derived normalized thresholds are genuinely useful when the oscillation amplitude of `eVarIx` changes substantially across parameter sweeps. In those workflows, absolute thresholds can become brittle or require retuning per regime.
- That utility is conceptually separate from Schmitt hysteresis. A warmup-derived normalized threshold can still be a single-boundary directional crossing, which is why `Observer.normalized_threshold_crossing` now exists as a separate lean family.
- The current `threshold_2` implementation combines three ideas at once: warmup-derived fractional placement, Schmitt-style hysteresis, and optional derivative gates.
- The derivative gates appear most justifiable on noisy or stochastic traces where value thresholds alone still overcount shallow or ripple-driven crossings. They are not the core semantic reason to keep `threshold_2`.
- This is the clearest evidence that the threshold family needed to be factored more cleanly rather than adding more options to `Observer.threshold_crossing` first.

### Constraint from the existing stub

- The current semantic `observer_threshold_crossing.clh` file is the lean proof target for absolute threshold observers, while the older `observer_threshold_1.clh` file is only historical evidence.
- It should not be treated as the implementation blueprint for the next slice.
- As written, it bakes in all-state trajectory extrema, event snapshots, event-window extrema, local-max and local-min side channels, and threshold-state machinery into one monolithic layout.
- Wiring that stub directly would move clODE away from the cleaner declaration boundary established by the summary slice rather than extending it.

### Why `neighbourhood_1` is a warning sign

- `nhood1` shows that “one pass” is not automatically the simplest design.
- Because its anchor point and normalization basis are discovered during the same pass that uses them for events, the trigger geometry can move as the solve proceeds.
- That is a real pitfall for any future declaration model: families with online-discovered geometry should probably not be treated as if they were just runtime-threshold variants of static one-pass triggers.
- Given the current weak workflow evidence, `nhood1` should stay under explicit keep-or-retire review rather than serving as the next abstraction target.

### Why two-pass families should not be grouped only by warmup

- `threshold_2` and `nhood2` both need a warmup pass, but they use that warmup for different reasons.
- `threshold_2` derives scalar trigger thresholds from extrema on one observed variable.
- `nhood2` derives normalization ranges from warmup, then discovers the anchor point `x0` online in the second pass only after a threshold-qualified crossing in `eVarIx`.
- A future declaration model should therefore separate warmup-derived geometry from second-pass anchor discovery instead of flattening both into one generic “two-pass” flag.

### Likely next design split

The current code suggests five observer-design axes that matter more than the current monolithic built-in names:

1. trigger geometry
2. parameter source: static, live-pass discovery, or warmup-derived
3. running-summary bundle
4. retained sparse-output policy
5. selected variable scope for any all-state or all-aux summary bundles

That does not mean every axis should become a user-facing object immediately. It does mean future observer-family declarations should probably be organized around those axes instead of around one opaque built-in selector.

The strongest immediate application of those axes is a deterministic static-trigger family with:

1. a scalar crossing trigger on `eVarIx` in absolute state-variable units
2. explicit direction or polarity selection
3. a separately selected `fVarIx` for event-conditioned summaries
4. independently chosen retained sparse-output policy
5. optional running-summary bundles that do not force a rebuild unless they change the schema

The minimum absolute-threshold slice is therefore a useful design checkpoint because it should need less build specialization than the summary family, not more.

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

- generalized scalar trigger functions are a good long-term umbrella, but the immediate implementation target should stay smaller: one-variable absolute threshold crossing with direction selection before any broader trigger-function DSL or code-generation surface is introduced.

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
