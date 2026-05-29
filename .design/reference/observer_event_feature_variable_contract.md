# Observer Event and Feature Variable Contract (eVarIx vs fVarIx)

Purpose: clarify the distinction between event-trigger variables and feature-measurement variables, and establish default resolution semantics.
Read when: configuring observer trigger and readout behavior, implementing event-detection or feature-extraction extensions, or understanding why observers can measure features in a variable different from their event trigger.
Update when: observer configuration semantics change, new observer families are introduced, or the default resolution strategy is revised.

## Overview

Observers expose two independent variable selectors:

- **eVarIx** ("event variable index"): The state variable on which to test for event conditions (e.g., threshold crossings, extrema)
- **fVarIx** ("feature variable index"): The state variable on which to measure features such as extrema timing, amplitude, and statistics

This separation enables powerful use cases while maintaining clarity about what each observer measures.

## Use Cases and Design Rationale

### 1. Independent Event Trigger and Feature Measurement

A user may wish to:
- Trigger events on one variable's oscillations (e.g., a fast feedback variable)
- Measure amplitudes and extrema in a different variable (e.g., a slower state that drives the behavior)

**Example**: In a coupled oscillator system with a slow variable $x_1$ and a fast variable $x_2$:
- Trigger local maxima events in $x_2$ (fast oscillator) to mark cycle boundaries
- Measure the amplitude of $x_1$ (slow variable) during each cycle

### 2. Simplifying Configuration Without Losing Power

**Threshold and Schmitt observers** benefit most from this flexibility:
- `threshold_crossing(eVarIx=0, fVarIx=0)`: Cross a threshold on $x_0$, measure amplitude in $x_0$
- `threshold_crossing(eVarIx=0, fVarIx=1)`: Cross a threshold on $x_0$, measure amplitude in $x_1$

Users can opt into advanced configurations without cluttering the default case.

### 3. Local Extrema Observers

**local_maximum** and **local_extremum** (now maxima-only) both use fVarIx to specify the variable on which to detect and measure extrema:
- `local_maximum(fVarIx=0)`: Trigger and measure in $x_0$
- `local_maximum(fVarIx=1)`: Trigger and measure in $x_1$

These observers do not currently use eVarIx (it is ignored), but the configuration surface is consistent.

### 4. Neighborhood Observers

**neighborhood_return** is unique: it uses all state variables to determine entry and exit from a neighborhood ball, so fVarIx specifies the variable used to detect the anchor point (downward threshold crossing) and is independent of the multi-variable neighborhood geometry.

## Default Resolution

### Recommended User-Facing Behavior

1. **If user specifies only one variable parameter** (current practice for simplicity):
   - Default both eVarIx and fVarIx to that variable index
   - This preserves backward compatibility and matches user intent in 95% of cases

2. **If user specifies event and feature variables separately** (future power-user case):
   - Respect both selectively
   - Documentation should clarify that some observers (e.g., neighborhood) may override one or both

### Implementation in Python Config Types

Observer config types should expose a single `variable` or `feature_variable` parameter by default:

```python
ThresholdCrossingConfig(event_var="x0", feature_var="x1")  # optional explicit separation
ThresholdCrossingConfig(variable="x0")                      # shorthand: both → "x0"
```

When converting to `ObserverRuntimeSettings`:
- Map `variable` → both `e_var_ix` and `f_var_ix` if both refer to the same variable
- Honor explicit separation if both are specified

## Family-Specific Semantics

### Threshold Crossing and Normalized Threshold Crossing

- **eVarIx**: Variable on which threshold crossings are detected
- **fVarIx**: Variable on which amplitude and other features are measured
- **Default**: Both point to the same variable (typical case)
- **Common**: `eVarIx = fVarIx` (cross a threshold on $x_0$, measure amplitude in $x_0$)

### Schmitt Trigger and Normalized Schmitt Trigger

- **eVarIx**: Variable on which hysteresis thresholds are tested (ignored if Schmitt-specific slope gates are used; recommended future path)
- **fVarIx**: Variable on which amplitude and extrema are measured
- **Default**: Both point to the same variable
- **Note**: Schmitt-specific readouts (up/down durations, duty cycle) are measured in $x_{eVarIx}$; trajectory and oscillation stats are measured in $x_{fVarIx}$

### Local Extremum and Local Maximum

- **eVarIx**: Ignored (not used for event detection)
- **fVarIx**: Variable on which local extrema are detected and measured
- **Default**: fVarIx = 0 (or user-specified)
- **Note**: These observers are extrema-first, not threshold-first; eVarIx is present for consistency but not actively used

### Neighborhood Return

- **eVarIx**: Ignored (anchor detected via downward threshold crossing in fVarIx)
- **fVarIx**: Variable on which the anchor threshold is tested; neighborhood geometry uses all variables
- **Default**: fVarIx = 0 (or user-specified)
- **Note**: The neighborhood ball is defined in normalized multi-variable space; fVarIx only specifies the anchor detection variable

## Migration and Compatibility

### Current Implementation

As of Phase 4 (local_extremum/local_maximum merge):
- All semantic event observers accept `f_var_ix` (feature/feature variable index) at the kernel level
- Python config types map `variable` → `f_var_ix` and `event_direction` appropriately
- eVarIx is currently always set to 0 or derived from compatibility surfaces (ObserverParams); it is not actively used in threshold or extremum families yet

### Future Direction (Phase 5+)

Once all semantic families have full parity:
- Introduce explicit `event_var` and `feature_var` config parameters on threshold and Schmitt families
- Expose eVarIx routing at the Python level through refined config classes
- Document this advanced configuration in power-user sections of the API docs
- Keep default behavior unchanged: single `variable` parameter maps to both eVarIx and fVarIx

### Backward Compatibility

- Existing code that specifies `threshold_crossing(variable="x0")` will continue to work
- Both eVarIx and fVarIx will be set to the same index
- No breaking changes to public API

## Recommendations for New Observers

When implementing a new observer family, consider:

1. **Does this observer need independent event and feature variables?**
   - Yes (e.g., threshold-based): design with separate eVarIx and fVarIx routing
   - No (e.g., neighborhood geometry): use fVarIx for all variable selection; document that eVarIx is ignored

2. **What is the semantic meaning?**
   - Event-first (threshold, Schmitt): eVarIx is primary; fVarIx is secondary for features
   - Extrema-first (local extrema): fVarIx is primary; eVarIx is not used
   - Geometry-first (neighborhood): fVarIx specifies anchor detection; all variables define geometry

3. **How should defaults work?**
   - Provide a single `variable` parameter in the Python config class
   - Let `to_observer_params()` set both indices to that variable for clarity
   - Document the advanced case where users can override one or both if needed

## Observed Patterns and Conventions

- All semantic event observers already track fVarIx consistently
- eVarIx is present in ObserverRuntimeSettings but under-utilized in current semantic families
- Threshold families are the natural candidates for future eVarIx/fVarIx separation
- Extrema families intentionally use fVarIx only; no event-trigger variable distinction is meaningful (extrema are always detected in the same variable they are measured)
