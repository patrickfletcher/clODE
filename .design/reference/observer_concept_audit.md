# Observer Concept Audit

Purpose: record the current, implementation-verified observer model in clODE and identify the explicit open design questions.
Read when: changing observer families, observer config surfaces, observer runtime settings, or observer output/readout semantics.
Update when: observer family behavior, observer config classes, kernel observer contracts, or verified observer tests change.

## Bottom line

- Observers are compile-time selected kernel families with persistent per-work-item observer state and a fixed readout schema resolved on the Python side.
- Canonical semantic observer families are `summary`, `threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, `local_max`, and `normalized_neighborhood_return`.
- Family-specific semantic config classes are live for threshold, Schmitt, local-maximum, and neighborhood-return observers.
- Compatibility surfaces (`ObserverParams`, legacy `observer_*` constructor kwargs, and `set_observer_parameters(...)`) remain supported.
- The active open design question is the shared seam for recurring event-trigger controls plus family-local readouts tracked in `.design/next_pr.md`.

## Fast path

- For current readout schema facts, use `.design/reference/observer_readout_audit.md`.
- For canonical versus compatibility policy, use `.design/reference/compatibility_boundary_audit.md`.
- For `e_var_ix` and `f_var_ix` routing, use `.design/reference/observer_event_feature_variable_contract.md`.
- For active implementation scope, use `.design/next_pr.md`.
- Historical rationale archived from this note is in `.design/archived/observer_concept_cleanup_2026_05_31/observer_concept_audit.md`.

## Verified implementation state

### Python observer surfaces

- `clode.observers.types.Observer` defines one canonical enum entry per semantic family.
- `ThresholdCrossingConfig` supports `threshold_crossing` and `normalized_threshold_crossing`.
- `SchmittTriggerConfig` supports `schmitt_trigger` and `normalized_schmitt_trigger`.
- `LocalMaximumConfig` supports `local_max` and enforces one-variable maxima detection (`e_var_ix == f_var_ix` through config mapping).
- `NeighborhoodReturnConfig` supports `normalized_neighborhood_return`.
- `FeatureSimulator.set_observer_configuration(...)` and `get_observer_configuration(...)` round-trip semantic config objects for the active observer family.
- `FeatureSimulator.set_observer_parameters(...)` remains the compatibility path for broad parameter updates.

### Build and runtime model

- Observer definitions resolve through `clode/observers/_definitions.py` into `ResolvedObserverSpec`.
- Readout feature names and observer-state layouts are resolved per observer definition (and summary variant when applicable).
- `SourceBuilder` compiles features kernels with one observer define plus resolved observer build variant/preamble.
- `OpenCLFeatureExecutor` runs one resolved observer spec per built features program.

### Kernel observer contract

`clode/kernels/observers.cl` defines the shared observer method lifecycle and includes family CLH implementations:

- `initializeObserverState(...)`
- `warmupObserverState(...)` (for two-pass families)
- `initializeEventDetector(...)`
- `updateObserverState(...)`
- `eventFunction(...)`
- `computeEventFeatures(...)`
- `finalizeFeatures(...)`

Shared accepted-step history helpers (`K=2` and `K=3`, including per-variable variants) are live in this file.

### Family behavior summary

| Family | Pass mode | Event detection channel | Notes |
| --- | --- | --- | --- |
| `Observer.summary` | one pass | none | summary reductions only |
| `Observer.threshold_crossing` | one pass | `eVarIx` | absolute threshold crossing |
| `Observer.normalized_threshold_crossing` | two pass | `eVarIx` | warmup-derived threshold |
| `Observer.schmitt_trigger` | one pass | `eVarIx` | absolute Schmitt transitions |
| `Observer.normalized_schmitt_trigger` | two pass | `eVarIx` | warmup-derived Schmitt transitions |
| `Observer.local_max` | one pass | `fVarIx` | local-extrema event streams, IMI/amplitude outputs |
| `Observer.normalized_neighborhood_return` | two pass | anchor/threshold on `eVarIx`; extrema tracking uses `fVarIx` | normalized neighborhood-exit detector |

## Verified tests and evidence

The current tree includes direct tests for:

- semantic config round-trips for threshold and Schmitt families, including normalized Schmitt
- Schmitt threshold validation (`x_up_threshold >= x_down_threshold`)
- rejection of derivative-threshold compatibility knobs on semantic Schmitt families
- normalized-threshold range validation and positive neighborhood radius validation
- feature-name/readout expectations for threshold and Schmitt families
- local-max three-sample timestamp refinement behavior

See:

- `test/test_simulation_contracts.py`
- `test/test_features.py`

## Active design work (explicitly not implementation facts)

- Bundle packaging remains open: should recurring event-trigger controls and family-local readouts be shared through one config seam or remain family-local?
- Schmitt duration outputs (`up duration`, `down duration`, `duty`) are now part of the current canonical Schmitt readout schemas; the remaining open question is whether any broader oscillation bundle seam should absorb them.
- `min_amp` and `max_event_count` now have clearer intended roles than the broader `ObserverParams` bundle suggests, and both now have live semantic homes across the current event-triggering families; the remaining open design question is packaging rather than adoption.

These are planning/implementation-tracking items, not all landed behavior yet; active tracking surfaces remain `.design/next_pr.md` and `.design/ideas.md`.

## Evidence anchors

- `clode/observers/types.py`
- `clode/observers/_definitions.py`
- `clode/simulation/features.py`
- `clode/_opencl/source_builder.py`
- `clode/_opencl/executors.py`
- `clode/kernels/observers.cl`
- `clode/kernels/observers/observer_threshold_crossing.clh`
- `clode/kernels/observers/observer_schmitt_trigger.clh`
- `clode/kernels/observers/observer_local_maximum.clh`
- `clode/kernels/observers/observer_normalized_neighborhood_return.clh`
- `test/test_simulation_contracts.py`
- `test/test_features.py`
