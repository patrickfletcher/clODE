# Observer Concept Audit

Purpose: audit what clODE observers currently mean in code and docs, compare that model to standard ODE-package event functions, and record what the landed summary-observer slice taught about declaration structure, build specialization, and runtime configuration.
Read when: planning observer refactors, comparing clODE to `solve_ivp`-style events, or deciding how Python should model event detection, running feature reduction, and observer outputs.
Update when: the observer concept changes materially, a new observer-family declaration slice lands, or trajectory/output policy grows close enough to observers that the boundary needs to be restated.

## Bottom line

- clODE observers are not just SciPy-style event functions. They are compile-time-selected, stateful, on-device feature pipelines that may detect events, accumulate online statistics, retain sparse event outputs, and emit a structured readout without storing full trajectories.
- The first declaration-model slice is now landed for summary-only reducers: `Observer.summary` plus `SummaryObserverSelection` resolve to a build-specialized summary variant, while `Observer.basic` and `Observer.basic_all_variables` remain compatibility presets over that family.
- That landed slice exposed a reusable pattern: a Python-side declaration resolves to one concrete spec that determines feature names, persistent layout, and build identity together, while scalar runtime knobs remain in the runtime-settings path.
- The first deterministic static-trigger slice is now also landed: `Observer.threshold_1` provides one-pass absolute threshold crossings with a lean timestamp-plus-count readout and runtime direction selection.
- The current `threshold_2` semantics are more specific than the name suggests: it is the live warmup-derived fractional threshold family, with Schmitt-style up/down state and optional slope gates layered on top.
- The threshold-family semantic config seam is now also landed for the current built-ins: threshold crossing and Schmitt triggering now have separate Python-side config objects with only the relevant knobs, while `ObserverParams` and the `observer_*` keyword arguments remain compatibility surfaces.
- The heavier event observers still bundle too many concerns under one built-in mode: trigger semantics, warmup or two-pass behavior, persistent state layout, feature-output schema, and retained event-output policy.
- That threshold slice confirmed the intended boundary from the opposite side of the summary family: scalar trigger controls such as crossing direction can stay runtime-side when the readout schema remains fixed.
- The next useful threshold proof target is therefore no longer the naming or parameter seam itself. That seam is now good enough to carry the next missing family: a lean warmup-derived fractional directional-threshold observer.
- `neighbourhood_1` currently looks more like a keep-or-retire case than a good template for future observer generalization unless a more deterministic workflow emerges.

## Current live model

### Python surface

- `clode.observers.types.Observer` still selects one built-in family or compatibility preset.
- The threshold families now also expose preferred semantic aliases at that surface: `Observer.threshold_crossing` for the current absolute one-boundary threshold observer and `Observer.schmitt_trigger` for the current warmup-derived hysteresis family. The older `threshold_1` and `threshold_2` names remain compatibility spellings.
- `ThresholdCrossingConfig` and `SchmittTriggerConfig` are now the first observer-family-specific semantic config objects beyond the summary selection surface.
- `Observer.summary` plus `SummaryObserverSelection` are now the first explicit family-specific declaration surface for observer feature selection.
- `ObserverParams` remains the public compatibility bundle, while `ObserverRuntimeSettings` and `EventOutputSettings` carry the narrower internal runtime settings and retained-event policy.
- The split between `e_var_ix` and `f_var_ix` is semantically important, not historical residue: users often need to trigger on one variable while measuring features on another, such as a slow calcium-like variable for event geometry and a faster voltage-like variable for amplitudes or peak counts.
- `ObserverDefinition` in `clode/observers/_definitions.py` still declares a build define, a two-pass flag, feature-name generation, and a matched persistent or event layout.
- The summary slice added a first family-specific resolved declaration path, but the definition layer still does not model event condition, running reducer family, variable-selection intent, or retained event-output intent as reusable concepts across the heavier event families.
- `FeatureSimulator` still routes most observer configuration through one observer name plus broad `observer_*` compatibility arguments, with `summary_selection` as the first narrower family-specific override.

### Runtime and build surface

- `SourceBuilder` still emits exactly one observer define for a features program build.
- `OpenCLFeatureExecutor` still owns one resolved observer spec, one observer-runtime-settings buffer, one observer-state buffer, and one feature buffer for that build.
- The landed summary slice extended `ResolvedObserverSpec` with a build variant and injected source preamble, and `BuildKey` now distinguishes observer variants beyond just the observer define plus `N_STORE_EVENTS`.
- `ObserverOutput` remains a string-keyed structured-array readout whose schema is fixed by the resolved observer spec.

### Kernel surface

- `clode/kernels/observers.cl` still expects one `ObserverState` type and one fixed function family: initialize, optional warmup, update, event test, event-feature computation, finalization, and continuation cleanup.
- The summary slice moved `basic` and `basicall` behind one `observer_summary.clh` template whose array extents and output schema are specialized by injected preamble macros.
- `localmax` and `nhood1` remain one-pass event detectors.
- `thresh2` and `nhood2` remain two-pass event detectors with a warmup solve.
- Heavy event observers still allocate state arrays scaled by `N_VAR`, `N_AUX`, and sometimes `N_STORE_EVENTS`.

## Current built-in observer matrix

Pass count alone is not enough to classify the current observers. The code currently exposes three parameter-source classes:

- no learned internal parameters at all
- on-the-fly internal-parameter discovery during the live pass
- warmup-derived internal parameters from a first pass

That distinction matters because it changes what can stay in runtime settings, what likely belongs in a family-specific declaration, and what makes a trigger geometry stable or unstable over the solve window.

| Public modes | Trigger class | Pass structure | Internal-parameter source | First warmup pass required? | Retained sparse output | Persistent summary scope today | Design notes |
| --- | --- | --- | --- | --- | --- | --- | --- |
| `Observer.summary`, `Observer.basic`, `Observer.basic_all_variables` | none; online reduction only | one pass | none beyond the resolved declaration and scalar runtime selectors | no | none | selected summary subset for `summary`; preset subset for `basic`; all-state plus all-aux preset for `basic_all_variables` | cleanest family; selection changes layout and schema directly |
| `Observer.local_max` | local maximum in `fVarIx` via derivative sign change with three-sample refinement | one pass | static runtime selectors only | no | local max/min timestamps and values up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle | clean deterministic event family; a mirrored local-min polarity likely belongs in the same future extremum family |
| `Observer.neighbourhood_1` | entry into an `nHoodRadius` ball around `x0` | one pass | `x0` discovered online from the first local minimum in `eVarIx`; normalization range keeps evolving during the same pass | no formal warmup, but yes to live-pass parameter discovery | none beyond event count | all-state plus all-aux summary bundle | awkward family: trigger geometry depends on mutable normalization data during the live pass; weak evidence for keeping it as a future template |
| `Observer.neighbourhood_2` | exit from an `nHoodRadius` ball around `x0` in normalized state space | two pass | warmup-derived trajectory ranges and threshold; `x0` is then pinned in the live pass as the first threshold-qualified point on `eVarIx` | yes | exit timestamps up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle, plus range and `x0` readout | better behaved than `nhood1`; useful as a specialized periodicity trigger for complex limit cycles because it requires return to a specific ball in full state space |
| `Observer.threshold_1` | absolute threshold crossing on `eVarIx` with runtime direction selection | one pass | static runtime selectors only | no | one threshold-event timestamp stream up to `N_STORE_EVENTS` | lean event stream only; no all-state summary bundle | confirms that a fixed-schema event family can keep threshold value and direction in runtime settings |
| `Observer.threshold_2` | warmup-derived fractional up/down thresholds on `eVarIx` with Schmitt-style up/down state; local maxima in `fVarIx` contribute derived summaries | two pass | warmup-derived fractions of global `x` and `dx` ranges on `eVarIx` | yes | up/down transition timestamps up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle, plus period, duty, and active-dip bundle | valuable when oscillation amplitudes change across a sweep; slope gates are mainly a noise-oriented extra filter rather than the defining semantics |

## Near-term candidate family ranking

| Candidate family | Current evidence | Near-term fit | Design read |
| --- | --- | --- | --- |
| Absolute threshold crossing in state-variable units | now landed as `Observer.threshold_1` with a lean event stream and runtime direction selection | landed first proof target | confirms that not every observer-family extension needs schema-changing specialization |
| Warmup-derived fractional threshold crossing | `threshold_2` already proves the value of warmup-relative thresholds across sweeps, while `threshold_1` proves lean directional threshold semantics | best next proof target | separates threshold parameterization from Schmitt topology instead of bundling both under one family |
| Local extremum family | `localmax` is live, and the shared three-sample extremum helpers already support symmetric max/min geometry | strong adjacent alternative | likely one polarity-selectable family rather than separate ad hoc `localmax` and `localmin` observers |
| Generalized hyperplane crossing | conceptually matches the useful part of `solve_ivp`-style scalar events | later | promising long-term trigger umbrella, but it needs a declaration for scalar trigger functions or hyperplane parameters that the current runtime path does not yet package cleanly |
| Neighborhood return triggers | `nhood1` and `nhood2` are live, and `nhood2` has a substantially better-behaved trigger because warmup pins the normalization before the live pass anchors `x0` below a threshold on `eVarIx` | specialized, not next | keep `nhood2` as a specialized periodicity detector for workflows where simple threshold crossings or local extrema are too permissive; `nhood1` remains the clearest retirement candidate |

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
- `threshold_1` currently occupies the absolute plus single-boundary corner.
- `threshold_2` currently occupies the warmup-derived plus Schmitt-style corner.
- The missing obvious corner is a warmup-derived fractional directional threshold crossing with the same lean readout style as `threshold_1`.
- An absolute Schmitt family is plausible later, but it is not the most informative next proof target because the repo already has one clean absolute threshold family and one warmup-derived hysteresis family.

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
- Once that semantic config seam exists, the missing warmup-derived directional-threshold family can land without worsening the current naming debt.

### What the landed threshold-family config seam proved

- Family-specific config objects are enough to improve the UX without changing the executor or kernel contracts. The Python side can adapt semantic configs into the existing normalized runtime settings.
- Some knobs can still recur across families, such as `min_amp`, without forcing one giant public bundle. The better rule is to surface a knob where it is active for that family, not where it might someday become reusable elsewhere.
- The broad `ObserverParams` bundle and `observer_*` keywords are still useful compatibility layers, but they no longer need to define the preferred semantic vocabulary for threshold families.

### Why `neighbourhood_2` remains valuable even if it is not next

- `nhood2` is materially better behaved than `nhood1` because the warmup pass fixes the trajectory ranges first, and the live pass only accepts `x0` after `x[eVarIx]` drops below a threshold derived from that warmup.
- That matters when `eVarIx` is a slow variable used to identify a meaningful phase of a bursting or otherwise multiscale oscillation.
- The family is therefore not just a threshold surrogate. Its intended semantics are to define a small neighborhood around a good full-state anchor point and then detect departure after a return to that neighborhood.
- That can be a more robust periodicity detector than simple threshold crossings or local-extremum triggers when complex limit cycles revisit similar scalar values at several distinct phases.
- The design conclusion is narrower than “neighborhood triggers are bad.” The current conclusion is that `nhood2` is a plausible specialized family, while `nhood1` is the weak link and should not drive the general observer abstraction.

## Design implications from the trigger matrix

### Why parameter source is a first-class design axis

- Pass count alone does not explain what kind of declaration a family needs.
- The more useful axis is where the trigger geometry gets its internal parameters: nowhere, from static runtime knobs, from live-pass discovery, or from a warmup pass.
- That axis should likely sit alongside trigger geometry, running-summary groups, and retained sparse-output policy in the next observer design.

### Why `eVarIx` versus `fVarIx` must stay explicit

- The observer catalog should not collapse event geometry and feature readout onto one implied primary variable.
- A real workflow is to trigger on a slower or more reliable reference variable while measuring period, amplitude, maxima, or other summaries on a different faster variable.
- That means future threshold, extremum, and hyperplane-trigger families should preserve separate trigger-variable and feature-variable roles, even if they default to the same variable in simpler use cases.

### Why absolute threshold crossing is the best next counterexample to summary-style specialization

- The summary family proved the shape-changing side of the declaration boundary: selected reductions change layout and schema, so they must be build-specialized.
- An absolute-threshold family is the opposite kind of test: threshold values, hysteresis values, and crossing direction are scalar behavior knobs that should usually stay in runtime settings because they do not inherently change observer-state layout.
- That makes a one-pass absolute-threshold family the best next proof target for checking whether the declaration model can stay narrow when the readout schema is fixed.
- It is also more deterministic and more intuitive to users than the current neighborhood-return ideas.

### What the landed `threshold_1` slice taught

- A lean event family can keep a fixed output schema while still exposing a meaningful new runtime control, here `EventDirection`.
- The first slice did not need build-signature changes beyond selecting the observer family itself: changing absolute threshold values or crossing direction does not rebuild the program.
- Rewriting the old `observer_threshold_1.clh` sketch down to one timestamp stream plus count was the right call. The monolithic all-state snapshot design was not needed to get useful behavior.
- The first landed slice still leaves one genuine follow-on question open: how much more event-conditioned readout can be added before the family stops being runtime-only and starts needing declaration-level schema control.

## What the minimum absolute-threshold family proved

The landed `threshold_1` slice did not clone `threshold_2`'s full readout bundle into a one-pass shell. That smaller proof target was the right call.

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

## What `threshold_2` still teaches that `threshold_1` does not

- Warmup-derived fractional thresholds are genuinely useful when the oscillation amplitude of `eVarIx` changes substantially across parameter sweeps. In those workflows, absolute thresholds can become brittle or require retuning per regime.
- That utility is conceptually separate from Schmitt hysteresis. A fractional threshold can still be a single-boundary directional crossing.
- The current `threshold_2` implementation combines three ideas at once: warmup-derived fractional placement, Schmitt-style hysteresis, and optional derivative gates.
- The derivative gates appear most justifiable on noisy or stochastic traces where value thresholds alone still overcount shallow or ripple-driven crossings. They are not the core semantic reason to keep `threshold_2`.
- This is the clearest evidence that the next threshold follow-on should factor the threshold family more cleanly rather than adding more options to `threshold_1` first.

### Constraint from the existing stub

- The current `observer_threshold_1.clh` file is useful as evidence that absolute-threshold observers were already contemplated.
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

- The first deterministic static-trigger slice is now landed, and the threshold re-audit now has a clearer answer than before.
- The immediate threshold-design question is no longer how to expose the semantic split. That Python-side config seam is now landed for the current threshold families.
- The next concrete implementation target should be the lean two-pass fractional directional-threshold family, and it should land through the new threshold-family config surface rather than through the old generic parameter bundle.
- `threshold_2` should be documented as the current warmup-derived Schmitt-style family rather than as the generic threshold baseline.
- Optional slope gates should stay framed as a noise-oriented or stochastic-oriented refinement rather than as the core threshold abstraction.
- `localmax` should be reconsidered as one half of a likely extremum family with a max/min polarity switch.
- `nhood2` should remain documented as a specialized periodicity detector for complex limit cycles even though it is not the next proof target.
- `nhood1` should remain under explicit keep-or-retire review and should not define the future observer abstraction unless a clearer deterministic workflow emerges.
- Generalized `solve_ivp`-style scalar trigger functions or hyperplane crossings should remain a later generalization after the landed scalar-crossing slice and the next adjacent follow-on clarify what reusable trigger packaging is still missing.
- Trajectory variable-subset storage should remain separate output-policy work and stays tracked in `.design/ideas.md`.

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
- `clode/kernels/observers/observer_threshold_2.clh`
- `clode/kernels/features.cl`
- `clode/kernels/initializeObserver.cl`
- `docs/feature_extraction.md`
- `scipy.integrate.solve_ivp` event documentation
