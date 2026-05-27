# Observer Concept Audit

Purpose: audit what clODE observers currently mean in code and docs, compare that model to standard ODE-package event functions, and record what the landed summary-observer slice taught about declaration structure, build specialization, and runtime configuration.
Read when: planning observer refactors, comparing clODE to `solve_ivp`-style events, or deciding how Python should model event detection, running feature reduction, and observer outputs.
Update when: the observer concept changes materially, a new observer-family declaration slice lands, or trajectory/output policy grows close enough to observers that the boundary needs to be restated.

## Bottom line

- clODE observers are not just SciPy-style event functions. They are compile-time-selected, stateful, on-device feature pipelines that may detect events, accumulate online statistics, retain sparse event outputs, and emit a structured readout without storing full trajectories.
- The first declaration-model slice is now landed for summary-only reducers: `Observer.summary` plus `SummaryObserverSelection` resolve to a build-specialized summary variant, while `Observer.basic` and `Observer.basic_all_variables` remain compatibility presets over that family.
- That landed slice exposed a reusable pattern: a Python-side declaration resolves to one concrete spec that determines feature names, persistent layout, and build identity together, while scalar runtime knobs remain in the runtime-settings path.
- The heavier event observers still bundle too many concerns under one built-in mode: trigger semantics, warmup or two-pass behavior, persistent state layout, feature-output schema, and retained event-output policy.
- The next useful work is therefore an observer-design audit, not immediate expansion of the event-observer catalog. The open question is no longer whether finer-grained declaration is useful, but where the build-specialized boundary should sit for one-pass and two-pass event observers.

## Current live model

### Python surface

- `clode.observers.types.Observer` still selects one built-in family or compatibility preset.
- `Observer.summary` plus `SummaryObserverSelection` are now the first explicit family-specific declaration surface for observer feature selection.
- `ObserverParams` remains the public compatibility bundle, while `ObserverRuntimeSettings` and `EventOutputSettings` carry the narrower internal runtime settings and retained-event policy.
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
| `Observer.local_max` | local maximum in `fVarIx` via derivative sign change with three-sample refinement | one pass | static runtime selectors only | no | local max/min timestamps and values up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle | cleanest current event family; likely best candidate for the next declaration split |
| `Observer.neighbourhood_1` | entry into an `nHoodRadius` ball around `x0` | one pass | `x0` discovered online from the first local minimum in `eVarIx`; normalization range keeps evolving during the same pass | no formal warmup, but yes to live-pass parameter discovery | none beyond event count | all-state plus all-aux summary bundle | awkward family: trigger geometry depends on mutable normalization data during the live pass |
| `Observer.neighbourhood_2` | exit from an `nHoodRadius` ball around `x0` in normalized state space | two pass | warmup-derived trajectory ranges and threshold; `x0` still discovered online in the second pass | yes | exit timestamps up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle, plus range and `x0` readout | mixes warmup-derived normalization with live-pass anchor discovery |
| `Observer.threshold_2` | upward threshold crossing on `eVarIx` with optional slope gate and hysteresis; local maxima in `fVarIx` contribute derived summaries | two pass | warmup-derived `xUp`, `xDown`, `dxUp`, and `dxDown` from global extrema on `eVarIx` | yes | up/down transition timestamps up to `N_STORE_EVENTS` | all-state plus all-aux summary bundle, plus period, duty, and active-dip bundle | strongest candidate for separating trigger policy, retained output policy, and running-summary groups |

## Trigger and parameter taxonomy

### 1. Static trigger families

- Summary reducers and `local_max` do not need a first pass to discover internal trigger parameters.
- Their runtime knobs are mostly scalar selectors, thresholds, or counters.
- These families are the easiest place to separate trigger geometry from retained outputs or running summary groups.

### 2. Online-discovery trigger families

- `neighbourhood_1` does not run a formal warmup pass, but it still discovers internal trigger parameters during the live pass.
- The anchor point `x0` is set from the first local minimum in `eVarIx`, and the normalization basis in `eventFunction(...)` uses trajectory ranges that are still changing over time.
- That makes `nhood1` a distinct design case rather than just another one-pass detector.

### 3. Warmup-derived trigger families

- `threshold_2` and `neighbourhood_2` both need a first pass to establish internal geometry before the main event pass starts.
- That is a good fit for a two-pass family concept, but the current code still mixes warmup-derived geometry with live-pass feature bundles and retained-output policy under one built-in name.

## Design implications from the trigger matrix

### Why parameter source is a first-class design axis

- Pass count alone does not explain what kind of declaration a family needs.
- The more useful axis is where the trigger geometry gets its internal parameters: nowhere, from static runtime knobs, from live-pass discovery, or from a warmup pass.
- That axis should likely sit alongside trigger geometry, running-summary groups, and retained sparse-output policy in the next observer design.

### Why `neighbourhood_1` is a warning sign

- `nhood1` shows that “one pass” is not automatically the simplest design.
- Because its anchor point and normalization basis are discovered during the same pass that uses them for events, the trigger geometry can move as the solve proceeds.
- That is a real pitfall for any future declaration model: families with online-discovered geometry should probably not be treated as if they were just runtime-threshold variants of static one-pass triggers.

### Why two-pass families should not be grouped only by warmup

- `threshold_2` and `nhood2` both need a warmup pass, but they use that warmup for different reasons.
- `threshold_2` derives scalar trigger thresholds from extrema on one observed variable.
- `nhood2` derives normalization ranges from warmup, then still discovers the anchor point `x0` online in the second pass.
- A future declaration model should therefore separate warmup-derived geometry from second-pass anchor discovery instead of flattening both into one generic “two-pass” flag.

### Likely next design split

The current code suggests five observer-design axes that matter more than the current monolithic built-in names:

1. trigger geometry
2. parameter source: static, live-pass discovery, or warmup-derived
3. running-summary bundle
4. retained sparse-output policy
5. selected variable scope for any all-state or all-aux summary bundles

That does not mean every axis should become a user-facing object immediately. It does mean future observer-family declarations should probably be organized around those axes instead of around one opaque built-in selector.

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

Semantics that should not become clODE's primary observer model:

- arbitrary Python callbacks inside the timestep loop
- dense output as the defining event abstraction
- treating all feature extraction as just event detection

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

- The next useful work is an observer-design audit, not immediate expansion of the event-observer catalog.
- That audit should document which observer decisions belong in the build-specialized declaration path versus the runtime-settings path.
- It should compare one-pass and two-pass event observers against the landed summary pattern before deciding whether the next implementation slice should be count-only event modes, selected running-summary groups for event observers, or a smaller family-specific declaration cleanup.
- It should keep the current built-in observer matrix up to date so planning discussions stay grounded in current trigger geometry and parameter-source classes rather than in broad labels like "one-pass" or "two-pass" alone.
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
