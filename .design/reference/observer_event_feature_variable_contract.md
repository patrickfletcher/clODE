# Observer Event/Feature Variable Contract

Purpose: document the current, implementation-verified meaning of event-variable and feature-variable selectors across observer families.
Read when: changing observer config dataclasses, runtime settings routing, or docs that describe event/feature variable behavior.
Update when: any observer family changes how `e_var_ix`/`f_var_ix` are set or consumed.

## Bottom line

- `e_var_ix` and `f_var_ix` are both live today, but each family uses them differently.
- Threshold and Schmitt families detect events from `e_var_ix`.
- `local_max` detects extrema from `f_var_ix` (with `LocalMaximumConfig` currently forcing `e_var_ix == f_var_ix`).
- `normalized_neighborhood_return` uses both: anchor/warmup logic from `e_var_ix`, maxima/amplitude tracking from `f_var_ix`.
- Threshold, Schmitt, and neighborhood semantic configs now expose both `event_var` and `feature_var` where those families split trigger geometry from extrema/amplitude measurement; `LocalMaximumConfig` remains the one-channel exception for now.

## Current implementation contract

### Threshold crossing families

- Public config: `ThresholdCrossingConfig(event_var=..., feature_var=..., threshold=..., direction=...)`
- Routed fields:
  - `event_var` -> `e_var_ix`
  - `feature_var` -> `f_var_ix`
- Kernel behavior (`observer_threshold_crossing.clh`, normalized variant):
  - event detection and threshold interpolation use `eVarIx`
  - maxima/amplitude tracking uses `fVarIx`
  - `min_amp` gating still follows the trigger channel (`eVarIx`)

### Schmitt families

- Public config: `SchmittTriggerConfig(event_var=..., feature_var=..., x_up_threshold=..., x_down_threshold=...)`
- Routed fields:
  - `event_var` -> `e_var_ix`
  - `feature_var` -> `f_var_ix`
- Kernel behavior (`observer_schmitt_trigger.clh`, normalized variant):
  - Schmitt state transitions and threshold logic use `eVarIx`
  - maxima/amplitude and active-dip tracking use `fVarIx`
- Validation behavior:
  - semantic Schmitt rejects non-default `dx_*` thresholds
  - `x_up_threshold >= x_down_threshold` is required
  - normalized Schmitt additionally requires thresholds in `[0, 1]`

### Local maximum family

- Public config: `LocalMaximumConfig(event_var=...)`
- Routed fields:
  - config `event_var` sets one index, then maps to both `e_var_ix` and `f_var_ix`
- Kernel behavior (`observer_local_maximum.clh`):
  - event detection and event-value extraction use `fVarIx`
  - this is effectively one-variable extrema detection by design

### Normalized neighborhood return family

- Public config: `NeighborhoodReturnConfig(event_var=..., feature_var=..., anchor_threshold=..., radius=...)`
- Routed fields:
  - `event_var` -> `e_var_ix`
  - `feature_var` -> `f_var_ix`
- Kernel behavior (`observer_normalized_neighborhood_return.clh`):
  - warmup anchor threshold and anchor latch logic use `eVarIx`
  - maxima/amplitude tracking path uses `fVarIx`
  - semantic config now exposes the same trigger/measurement split as the live kernel

## Compatibility surface behavior

`ObserverParams`, constructor `observer_*` args, and `set_observer_parameters(...)` can still set both `event_var` and `feature_var` directly for all families. That path remains the broad compatibility adapter.

## Evidence anchors

- `clode/observers/types.py`
- `clode/simulation/features.py`
- `clode/kernels/observers/observer_threshold_crossing.clh`
- `clode/kernels/observers/observer_schmitt_trigger.clh`
- `clode/kernels/observers/observer_local_maximum.clh`
- `clode/kernels/observers/observer_normalized_neighborhood_return.clh`
- `test/test_simulation_contracts.py`
- `test/test_features.py`

## Decided direction and follow-through

Direction is now explicit:

- where kernels consume both `eVarIx` and `fVarIx` semantics, semantic config surfaces should expose both `event_var` and `feature_var`.

Current live split:

- `ThresholdCrossingConfig` now exposes `feature_var`; threshold event geometry and `min_amp` stay on `eVarIx`, while extrema/amplitude tracking uses `fVarIx`.
- `NeighborhoodReturnConfig` now exposes `feature_var`; `observer_normalized_neighborhood_return.clh` already consumes `fVarIx` for maxima/amplitude tracking.
