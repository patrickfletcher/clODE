# Observer Concept Audit

Purpose: audit what clODE observers currently mean in code and docs, compare that model to standard ODE-package event functions, and define the most plausible implementation path for finer-grained feature selection and memory-footprint control.
Read when: planning observer refactors, comparing clODE to `solve_ivp`-style events, or deciding how Python should model event detection, running feature reduction, and observer outputs.
Update when: the observer concept changes materially, a first observer-authoring slice lands, or trajectory/output policy grows close enough to observers that the boundary needs to be restated.

## Bottom line

- clODE observers are not just SciPy-style event functions. They are compile-time-selected, stateful, on-device feature pipelines that may detect events, accumulate online statistics, retain sparse event outputs, and emit a structured readout without storing full trajectories.
- The current implementation works, but one built-in observer name still selects too many concerns at once: trigger semantics, warmup or two-pass behavior, persistent state layout, feature-output schema, and retained event-output policy.
- That bundling makes extension hard and gives users coarse control over memory footprint. Aside from `basic` plus the generic `event_var` and `feature_var` selectors, most built-ins either track one distinguished variable or eagerly allocate and read out all state variables and auxiliary variables.
- The next useful implementation path is not a public custom-observer DSL. It is an internal Python-side declaration model that separates event condition, running reducers, variable selection, and output-retention policy while still compiling one supported observer template per feature build.
- User-selectable trajectory variable subsets are related output-policy work, but they should stay separate from this observer audit and remain on the trajectory backlog rather than being folded into the observer PR.

## Current live model

### Python surface

- `clode.observers.types.Observer` selects one built-in mode.
- `ObserverParams` remains the public compatibility bundle, while `ObserverRuntimeSettings` and `EventOutputSettings` carry the narrower internal runtime settings and retained-event policy.
- `ObserverDefinition` in `clode/observers/_definitions.py` declares a build define, a two-pass flag, feature-name generation, and a matched persistent or event layout.
- That definition layer does not yet model event condition, running-reducer family, variable-selection intent, or output-retention intent as separate concepts.
- `FeatureSimulator` still routes observer configuration through one observer name plus broad `observer_*` compatibility arguments.

### Runtime and build surface

- `SourceBuilder` emits exactly one observer define for a features program build.
- `OpenCLFeatureExecutor` owns one resolved observer spec, one observer-runtime-settings buffer, one observer-state buffer, and one feature buffer for that build.
- Changing the observer build define or `N_STORE_EVENTS` changes the build key.
- `ObserverOutput` is a string-keyed structured array readout whose schema is fixed by the selected observer definition.

### Kernel surface

- `clode/kernels/observers.cl` expects one `ObserverState` type and one fixed function family: initialize, optional warmup, update, event test, event-feature computation, finalization, and continuation cleanup.
- `basic` and `basicall` are summary reducers with no event trigger.
- `localmax` and `nhood1` are one-pass event detectors.
- `thresh2` and `nhood2` are two-pass event detectors with a warmup solve.
- Heavy event observers allocate state arrays scaled by `N_VAR`, `N_AUX`, and sometimes `N_STORE_EVENTS`.

## Relation to `solve_ivp` events

`scipy.integrate.solve_ivp` events are scalar root functions of `(t, y)` with optional `direction` and `terminal` semantics. The solver finds zero crossings, returns `t_events` and `y_events`, and otherwise leaves summary-feature logic to user code layered on top of the solver.

That model overlaps with clODE only at the event-trigger layer:

- `threshold_2`, `local_max`, and the neighborhood observers have an event-condition component that is broadly comparable to `solve_ivp` events.
- `basic` and `basicall` do not correspond to event functions at all; they are online reducers over the whole solve window.
- clODE observers additionally own persistent per-trajectory state, online summary updates, event-count or timestamp retention, and a final structured readout. Those are central to clODE's large-ensemble workflow and should not be reduced to plain callback-style events.

Useful semantics to borrow from `solve_ivp` as naming or modeling guidance:

- explicit distinction between event condition and event readout
- explicit directionality for crossings
- explicit terminal vs non-terminal event behavior
- clear statement that multiple crossings inside one coarse step may still be missed unless the kernel tracks more geometry

Semantics that should not become clODE's primary observer model:

- arbitrary Python callbacks inside the timestep loop
- dense output as the defining event abstraction
- treating all feature extraction as just event detection

## Current friction points

### 1. Observer identity is monolithic

- One observer name currently chooses trigger semantics, persistent state layout, readout schema, warmup policy, and build-time kernel selection all at once.
- That is manageable for a small built-in catalog, but it blocks finer-grained feature selection and makes built-in extension feel heavier than it needs to be.

### 2. Runtime settings are broad and partially inactive

- `ObserverRuntimeSettings` exposes one shared set of indices, thresholds, and counters for all built-ins.
- Many fields are irrelevant for a given observer, so the Python surface does not tell users which knobs are structural and which are inert for the current mode.

### 3. Feature-selection granularity is coarse

- `basic` tracks one selected feature variable.
- `basicall` tracks all state variables plus all auxiliary variables.
- The heavier event observers also carry all-state summary arrays, even when a user may only care about period, count, or extrema on one or two variables.
- There is no current way to ask for “max/min/mean of variables x and z only” or “count local extrema without storing the full all-variable summary bundle” through the Python definition layer.

### 4. Memory footprint follows observer mode more than requested readout

- `basicall`, `localmax`, `nhood1`, `nhood2`, and `thresh2` size persistent arrays directly from `N_VAR` and often `N_AUX`.
- Event-capable observers also scale state or retained outputs with `N_STORE_EVENTS`.
- Today the main user control over observer footprint is indirect: pick a lighter observer, reduce `max_event_timestamps`, or avoid auxiliary-heavy models.

### 5. Output schema is fixed and stringly typed

- `ObserverOutput` is useful for inspection, but the schema is still just an ordered list of strings paired with one flat structured array.
- That makes it easy to expose results, but harder to declare feature families, optional groups, or memory-affecting output policy on the Python side before the build happens.

### 6. Compile-time builds limit free-form composition

- The current build model is still a strength, not a bug: one compiled observer template per build keeps OpenCL layouts explicit and testable.
- It also means clODE should not promise arbitrary mix-and-match observer composition immediately. Any new user-facing selection surface needs to map to a small set of supported kernel templates and state layouts.

## Questions this audit should answer before implementation

- Which parts of the observer concept deserve separate Python-side declaration objects?
- Which current built-ins are really the same template family with different readout policies?
- Which user controls would actually reduce state size or output size, and which would only rename fields without changing footprint?
- How much `solve_ivp`-style event vocabulary should clODE adopt for triggers, directions, and terminal behavior without collapsing the broader observer concept into callback-style events?
- Which observer selections belong in the same workstream as trajectory-output policy, and which should remain separate?

## Recommended implementation path

### 1. Keep one compiled observer template per feature build in the near term

- Preserve the current compile-time observer selection model while the declaration layer is cleaned up.
- Treat the observer build choice as a template family, not yet as a free-form composition engine.

### 2. Split the internal declaration model into smaller concepts

The next Python-side observer model should make these concerns explicit even if they remain internal dataclasses at first:

- event condition or detection family
- running reducer or summary family
- variable-selection policy
- retained event-output policy
- warmup or two-pass requirement

This should happen in Python first, with `_opencl` consuming the declarations rather than inventing parallel structure.

### 3. Reclassify the current built-ins into a few template families

The existing catalog already suggests three coarse families:

- summary-only reducers
- one-pass event detectors
- two-pass event detectors

That family split is likely a better implementation anchor than keeping every current built-in as a totally separate conceptual island.

### 4. Add user-facing selection only where it maps to smaller supported templates

Likely high-value early controls include:

- selected variables for summary reducers rather than only one variable or all variables
- count-only vs timestamp-retaining event modes where the kernel template can truly drop state
- clearer distinction between event variable, feature variable, and any future selected summary-variable subset

Likely poor early targets include:

- arbitrary user-defined algebra over state variables
- unconstrained composition of every reducer with every event trigger
- public custom observer authoring before the built-in declaration model stabilizes

### 5. Use the `solve_ivp` comparison to sharpen terminology, not to copy the full model

- Treat `solve_ivp` events as the closest analogue for clODE event-condition semantics.
- Treat clODE observers as the larger envelope: event condition plus online reduction plus retained sparse output plus final readout.
- Keep trajectory output policy separate from that envelope.

## Recommended first proof target

The best first implementation slice is a unified summary-observer family behind the current `basic` and `basicall` surface.

Why this slice is the strongest:

- it is one-pass and event-free, so it exercises observer declaration and feature selection without reopening warmup passes, crossing geometry, or timestamp retention
- it targets the most obvious missing user control today: “track max/min/mean or slope summaries for these variables, not all variables”
- it is the clearest place where smaller selections can and should reduce both feature schema size and persistent observer-state size
- it preserves the heavier event observers while giving them a better declaration model to grow on later

What this slice should prove:

- `basic` and `basicall` are really one summary template family with different selection presets
- Python-side selection policy should be explicit for tracked state variables, tracked auxiliary variables, and whether slope extrema are included
- current public modes can remain compatibility presets while the internal model grows more expressive
- subset-driven state reduction likely needs family-specific configuration that is more precise than the current shared `ObserverRuntimeSettings` bundle

What this slice should not try to do yet:

- general event-detector composition
- public custom observer authoring
- trajectory variable-subset storage
- a new public observer DSL

## Consequences for the next implementation PR

- The first implementation slice should be a summary-family declaration cleanup, not a broad kernel rewrite.
- The validation target should be unifying `basic` and `basicall` under one summary-family declaration with a selected-variable policy and compatibility presets for the current public modes.
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
- `clode/kernels/observers.cl`
- `clode/kernels/features.cl`
- `clode/kernels/initializeObserver.cl`
- `docs/feature_extraction.md`
- `scipy.integrate.solve_ivp` event documentation
