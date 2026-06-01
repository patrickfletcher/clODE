# Observer Readout Audit

This document analyzes the readout families across legacy and semantic observers to establish architectural seams for selectable readout subsets.

## Fast path

- Treat `local_maximum` semantics as the canonical oscillation contract (`period`, `amplitude`, optional max/min event streams).
- Avoid adding new `local_max`-specific readout complexity while this contract is being standardized.
- Treat `local_max` as the public extrema family; keep compatibility-path discussion separate from core readout semantics.
- Keep helper-hardening work (compensated time/mean, interpolation variants) as non-blocking backlog.

## Readout Categorization

### Category 1: Summary State Readouts (Always Available)

These readouts apply to **all** observers and should be available uniformly:

- **Trajectory extrema**: `max {var}`, `min {var}` for each state variable
- **Trajectory extrema of derivatives**: `max d{var}/dt`, `min d{var}/dt` for each state variable
- **Trajectory mean (time-integrated)**: `mean {var}` for each state variable
- **Auxiliary trajectory extrema**: `max {aux}`, `min {aux}` for each auxiliary variable
- **Auxiliary trajectory mean**: `mean {aux}` for each auxiliary variable

These are independent of the event-detection mechanism and reflect the full trajectory statistics during observation. They are managed by helper functions already (`xTrajectoryMax`, `xTrajectoryMin`, etc.).

**Strategic Notes**:
- It is possible to compute max/min/mean period and amplitude without event functions by tracking successive local maxima/minima in a monitored variable.
- Canonical oscillation semantics: `period = tThisMax - tLastMax`; `amplitude = xThisMax - xLastMin`.

### Category 2: Event-Stream Readouts (Any Observer with Event Function)

These readouts are available for **any** observer that detects events (has non-null `eventFunction`):

- **Event times**: `event time 0`, `event time 1`, etc. (up to `N_STORE_EVENTS`)
- **Event count**: `event count`
- **Event state values** (future): `event value 0`, `event value 1`, etc. (not yet implemented, but intended for `solve_ivp` compatibility)

**Current observers with events**:
- `threshold_crossing` (1-pass, absolute)
- `normalized_threshold_crossing` (2-pass, warmup-derived)
- `schmitt_trigger` (1-pass, absolute)
- `normalized_schmitt_trigger` (2-pass, warmup-derived)
- `local_max` (1-pass)
- `normalized_neighborhood_return` (2-pass)

**Current storage**: Timestamps only; event state/value extraction is not yet standardized across families.

**Strategic Notes**:
- Standardize event-stream names across observers as `event time {event_ix}` and `event value {var_ix} {event_ix}`.

### Category 3: Period Statistics (Event-Detector Families Only — Universal Name, Context-Dependent Meaning)

Available for observers that define "periods" (inter-event intervals) between events:

- **Period statistics**: `max period`, `min period`, `mean period`

The meaning of "period" varies by event trigger geometry, but the name and computation are uniform:

| Observer | Period Meaning | Computation |
|----------|---------------|------------|
| `local_max` | Inter-extremum interval (time between consecutive maxima or minima) | Elapsed time between events |
| `threshold_crossing` | Inter-crossing interval (time between consecutive threshold crossings) | Elapsed time between events |
| `normalized_threshold_crossing` | Inter-crossing interval (time between consecutive crossings) | Elapsed time between events |
| `schmitt_trigger` | Inter-transition interval (time between consecutive up-transitions, typically) | Elapsed time between events |
| `normalized_schmitt_trigger` | Inter-transition interval (time between consecutive up-transitions) | Elapsed time between events |
| `normalized_neighborhood_return` | Inter-exit interval (time between consecutive neighborhood exits) | Elapsed time between events |

**Strategic Notes**:
- "Period" = "inter-event interval" universally; exact meaning depends on what the event trigger detects
- Computed using **elapsed time** (compensated for adaptive stepping), not wallclock time
- Some families compute period only when `eventcount > 1` (at least 2 events detected); single-event scenarios return ZERO
- Period is fundamentally meaningful for all non-summary event detectors and will be computed universally in Phase 2 onwards

### Category 4: Amplitude and Maxima Counting (Oscillation Readouts Across Event Observers)

#### Amplitude Statistics

Available for any observer with events, by tracking local maxima and minima in the monitored feature variable (`fVarIx`):

- **Amplitude statistics**: `max amplitude`, `min amplitude`, `mean amplitude`
  - Defined as: `local_max_value - last_local_min_value`
  - Measured in one variable specified by `fVarIx`

**Current observers with this**:
- `local_max` - computes amplitude

**Strategic Direction**: Amplitude should be standardized across event-observer families. Event detectors already maintain step history for period and maxima counting, so tracking opposite-polarity extrema is feasible and does not require coupling amplitude semantics to event polarity. This includes threshold, Schmitt, neighborhood, and local-extremum families.

#### Maxima Counting (Per-Event Local Maxima Count)

Available for observers where it is meaningful: **all event detectors EXCEPT local-extremum**.

- **Maxima count**: `max peak count`, `min peak count`, `mean peak count`
  - Counts local maxima between consecutive events
  - Measured in one variable specified by `fVarIx`

**Strategic Rationale**: We omit maxima counting from `local_max` to keep one canonical oscillation contract centered on `local_maximum` and avoid duplicate near-equivalent surfaces. Threshold, Schmitt, and neighborhood triggers still benefit from maxima counts because those event geometries are not already an extremum-stream abstraction.

**Current observers with this**:
- `threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, and `normalized_neighborhood_return` compute peak-count statistics where applicable.

**Planned (Phase 2+)**:
- Threshold and Schmitt semantic families will compute peak counts in Phase 2
- Neighborhood families (both semantic and legacy) may compute peak counts (Phase 2+, separate validation)

**Intentional Omission**:
- `local_max` does not compute peak-count statistics.

### Category 5: Schmitt-Specific Readouts (Schmitt Trigger Families Only)

Available only for Schmitt-style observers:

- **Up-state duration**: `max upDuration`, `min upDuration`, `mean upDuration`
  - Time spent in the "up" state before down-transition
- **Down-state duration**: `max downDuration`, `min downDuration`, `mean downDuration`
  - Time spent in the "down" state before up-transition
- **Duty cycle**: `max duty`, `min duty`, `mean duty`
  - `upDuration / period`

**Current observers**:
- `schmitt_trigger` and `normalized_schmitt_trigger` currently expose up/down transition timing streams and event count.

**Note**: These are heavily tied to the Schmitt state machine and are not applicable to threshold-crossing or extremum families.

### Category 6: Schmitt Legacy Extras (normalized_schmitt_trigger Only)

- **Active dip**: `max activeDip`, `min activeDip`, `mean activeDip`
  - Feature-variable value during down-state minus last recorded feature-variable minimum
  - Highly specific to the legacy workflow; **not recommended for standardization**

---

## Common Readout Processing Patterns

### Pattern 1: Running Min/Max Tracking

Used for computing extrema across events:
- Accumulate `max`, `min` for each event period
- Compute running mean after sufficient events

**Shared helper**: `runningMean()` (already available)

### Pattern 2: Three-Sample Refinement for Event Times

Used to refine event detection from sampled to interpolated times:
- `localMaximumFromThreeSamples()`
- `localMinimumFromThreeSamples()`

Already in shared helpers; used by `local_max`, `normalized_neighborhood_return`, and threshold/Schmitt families.

### Pattern 3: Time-Integrated Trajectory Statistics

Used for computing trajectory mean across the full time span:
- `runningMeanTime()` helper already in use

### Pattern 4: Event Stopping Control

Used to limit integration based on event count:
- `if (od->eventcount == op->maxEventCount) return true;`

This is a shared semantic: "stop after N events," analogous to `max_steps` and `max_store`.

---

## Proposed Helper Organization

### New/Consolidated Header: `clODE_observer_helpers.clh`

Consolidate and formalize shared readout computation into dedicated helpers. These should use best available:
- **Compensated elapsed time** (not wallclock `t`) for period and mean tracking
- **Bounded three-sample quadratic refinement** for event-time and extremum-value accuracy
- **`runningMeanTime()` helper** for time-integrated means with compensated time

```c
// Trajectory extent helpers (already exist inline; formalize for clarity)
inline void updateTrajectoryMaxMin(
    realtype *xMax, realtype *xMin,
    const realtype x_current,
    int variable_index
);

// Time-integrated mean: updates running mean using elapsed time (compensated)
inline void updateTrajectoryMean(
    realtype *xMean,
    const realtype x_current, const realtype stepDt,
    const realtype elapsedTotal
);

// Period tracking: accumulate inter-event intervals
inline void updatePeriodStats(
    realtype period[3],  // max/min/mean
    const realtype elapsedThisEvent,
    const realtype elapsedLastEvent,
    const int eventcount
);

// Amplitude tracking: store and compare extrema
inline void updateAmplitudeStats(
    realtype amplitude[3],  // max/min/mean
    const realtype thisExtremaValue,
    const realtype lastOppositeExtremaValue,
    const int eventcount
);

// Maxima counting: count local maxima between events (for non-extremum triggers)
inline void updateMaximaCountStats(
    realtype nMaxima[3],  // max/min/mean
    const unsigned int thisNMaxima,
    const int eventcount
);

// Auxiliary trajectory tracking (similar to state variable tracking)
inline void updateAuxTrajectoryStats(
    realtype *auxMax, realtype *auxMin, realtype *auxMean,
    const realtype aux_current, const realtype stepDt,
    const realtype elapsedTotal,
    int auxiliary_index
);
```

**Benefit**: Standardizes computation across families, reduces code duplication, and makes intent clear.

### Use of Existing Helpers

Leverage already-landed helpers:
- `advanceAcceptedStepHistory2()`, `advanceAcceptedStepHistory3()` and `ByVariable` variants for buffer management
- `localMaximumFromThreeSamples()`, `localMinimumFromThreeSamples()` for event refinement
- `runningMean()`, `runningMeanTime()` for aggregation
- Compensated-time tracking mechanisms already in use by adaptive steppers

---

## Recommended Readout Naming Convention

Use short, consistent names across families. Meaning is disambiguated by **context (which observer is selected)** and **documentation** rather than verbose naming.

### Standardized Names

| Category | Feature Names | Notes |
|----------|--------------|-------|
| Summary State | `max {var}`, `min {var}`, `mean {var}` | Available for all variables and aux |
| Summary Slope | `max d{var}/dt`, `min d{var}/dt` | Available for all variables |
| Event Streams | `event time 0`, `event time 1`, ... | Up to `N_STORE_EVENTS` |
| Event Count | `event count` | Total events detected |
| Period | `max period`, `min period`, `mean period` | Meaning documented per family |
| Amplitude | `max amplitude`, `min amplitude`, `mean amplitude` | Measured in `fVarIx` only |
| Peak Count | `max peak count`, `min peak count`, `mean peak count` | Local maxima between events |
| Up Duration | `max up duration`, `min up duration`, `mean up duration` | Schmitt-only |
| Down Duration | `max down duration`, `min down duration`, `mean down duration` | Schmitt-only |
| Duty Cycle | `max duty`, `min duty`, `mean duty` | Schmitt-only |

**Avoid in public API**:
- `IMI` (use "period" and document as "inter-maxima interval" for `local_maximum`)
- `activeDip` (too niche for legacy `normalized_schmitt_trigger`)
- `nhood center` (internal state; use structured variable selection instead)

---

## Current Observer Readout Inventory

**Important Context**: Semantic families (`threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, `local_max`, `normalized_neighborhood_return`) are **intentionally lean shells** designed to have feature readouts added consistently during planned phases. The gaps documented below are strategic, not oversights.

### `local_maximum` (Legacy, 1-Pass Extremum Detection)

**Available Readouts**:
- Summary trajectory stats: `max/min/mean {var}`, `max/min d{var}/dt`, `max/min/mean {aux}`
- Period statistics: `max/min/mean period` (IMI—inter-maxima interval)
- Amplitude: `max/min/mean amplitude` (max value - last min value in `fVarIx`)
- Maxima count: `max/min/mean peak count`
- Event streams: `event time 0-N`, event count (stored as `tMaxList`, `tMinList`, `xMaxList`, `xMinList`)

**Layout Complexity**: High (14 persistent fields + event arrays)

---

### `normalized_schmitt_trigger` (Legacy, 2-Pass Normalized Schmitt + Extras)

**Available Readouts**:
- Summary trajectory stats: `max/min/mean {var}`, `max/min d{var}/dt`, `max/min/mean {aux}`
- Period statistics: `max/min/mean period` (time between up-transitions)
- Up duration: `max/min/mean upDuration`
- Down duration: `max/min/mean downDuration`
- Duty cycle: `max/min/mean duty`
- Active dip: `max/min/mean activeDip` (feature-variable dep-state value - last min)
- Peak count: `max/min/mean peak count` (maxima per period in `fVarIx`)
- Event streams: up/down transition times (stored as `tUpTransition`, `tDownTransition`), event count

**Layout Complexity**: Very high (41+ persistent fields + event arrays)

---

### `normalized_neighborhood_return` (Legacy, 2-Pass Neighborhood Return)

**Available Readouts**:
- Summary trajectory stats: `max/min/mean {var}`, `min/max/range {var}`, `max/min d{var}/dt`, `max/min/mean {aux}`
  - Note: includes `range` (max - min) and `nhood center` (x0) for each variable
- Period statistics: `max/min/mean period` (time between neighborhood exits)
- Peak count: `max/min/mean peak count` (maxima per period in `fVarIx`)
- Event streams: neighborhood-exit times, event count

**Layout Complexity**: High (28+ persistent fields + event arrays)

---

### Semantic Threshold Family (`threshold_crossing`, `normalized_threshold_crossing`)

**Currently Available Readouts**:
- Event streams: event times, event count (minimal)

**Planned (Phase 2+)**:
- Summary trajectory stats: `max/min/mean {var}`, `max/min d{var}/dt`, `max/min/mean {aux}` (universal addition)
- Period statistics: `max/min/mean period` (inter-crossing intervals)
- Maxima count: `max/min/mean peak count` (local maxima between crossings)
- **Planned (Phase 3)**: Amplitude statistics `max/min/mean amplitude`

**Intentional Omissions** (by design as lean shell):
- Schmitt-specific readouts (not applicable to threshold crossing)

**Layout Complexity**: Very low currently (6 persistent fields + event array); will expand modestly in Phase 2

---

### Semantic Schmitt Family (`schmitt_trigger`, `normalized_schmitt_trigger`)

**Currently Available Readouts**:
- Event streams: up/down transition times, event count

**Planned (Phase 2+)**:
- Summary trajectory stats: `max/min/mean {var}`, `max/min d{var}/dt`, `max/min/mean {aux}` (universal addition)
- Period statistics: `max/min/mean period` (inter-transition intervals)
- Maxima count: `max/min/mean peak count` (local maxima between transitions)
- **Planned (Phase 3)**: Amplitude statistics `max/min/mean amplitude`
- **Future (Phase 4+)**: Up/down duration and duty cycle statistics (if feature parity with legacy `normalized_schmitt_trigger` is desired)

**Intentional Omissions** (by design as lean shell):
- None beyond Schmitt-specific extras deferred to later phase

**Layout Complexity**: Low currently (9 persistent fields + event arrays); will expand in Phase 2, optionally further in Phase 3-4

---

### Semantic Extremum Families (`local_max`, `normalized_neighborhood_return`)

**Currently Available Readouts**:
- Event streams: event times (and values for `local_max`), event count

**Planned (Phase 2)**:
- Summary trajectory stats: `max/min/mean {var}`, `max/min d{var}/dt`, `max/min/mean {aux}` (universal addition)
- Period statistics: `max/min/mean period` (inter-event intervals)

**Planned (Phase 3)**:
- `local_max`: Amplitude `max/min/mean amplitude` (local extremum value - last opposite-polarity extremum value)
- `normalized_neighborhood_return`: Amplitude `max/min/mean amplitude` (local oscillation amplitude in `fVarIx`, independent of trigger geometry)
- `local_max` still **omits maxima counting** (strategic decision: alternating extrema makes count uninformative)

**Intentional Omissions** (by design):
- Schmitt-specific readouts (not applicable to extremum or neighborhood triggers)
- Maxima counting for `local_max` (meaningless given trigger geometry)

**Layout Complexity**: Very low currently (4-8 persistent fields + event arrays); will expand modestly in Phase 2-3

---

## Design Findings and Strategic Decisions

**Key Strategic Context**: Semantic families were intentionally implemented as lean shells to have feature readouts built in *consistently once the basic event detection logic was established*. The goal is strategic planning based on theoretical possibilities, not on what is currently (sparsely) implemented.

### 1. Summary Trajectory Statistics Should Be Universal

**Finding**: Legacy families compute trajectory extrema; semantic families do not (by design—they are lean shells).

**Strategic Decision**: All event-detector observers should track trajectory statistics uniformly. This adds modest state overhead (`5*N_VAR + 3*N_AUX` floats) and ensures consistency.

**Implementation Approach**: Add trajectory-stats fields to semantic families in Phase 2. Use best available compensated-time and mean helpers.

### 2. Period (Inter-Event Interval) Should Be Universal

**Strategic Decision**: "Period" = "inter-event interval" across all event-detector families, with meaning varying by trigger geometry:

- **Local extremum**: Time between consecutive extrema (inter-maxima interval for maxima polarity, inter-minima for minima)
- **Threshold crossing**: Time between consecutive threshold crossings
- **Schmitt trigger**: Time between consecutive up-transitions (or down-transitions; convention is needed)
- **Neighborhood return**: Time between consecutive neighborhood-exit events

**Rationale**: This unifies the naming and clarifies that computation logic is similar (accumulated elapsed time between events) even though the trigger geometry differs.

**Implementation Approach**: Compute `period[3]` (max/min/mean) for all event-detector families in Phase 2. Use `elapsedTime` not wallclock `t` to account for time dilation in adaptive-step solvers.

### 3. Maxima Counting Is Only Meaningful for Non-Extremum Triggers

**Strategic Finding**: Maxima counting is most useful where event geometry is not itself a direct extremum stream.

**Strategic Decision**: Compute `nMaxima[3]` (max/min/mean count of local maxima per event) for threshold, Schmitt, and neighborhood families. Keep it out of `local_max` while the local_maximum-first simplification path is active.

**Implementation Approach**: Add maxima-counting logic to non-extremum families where triggered events are not themselves extrema detections.

### 4. Amplitude Should Be Universal Across Event Observers

**Strategic Finding**: Amplitude is defined as `local_max_value - last_local_min_value` and requires tracking both extrema polarities. That tracking can be maintained independently of event-trigger geometry.

**Strategic Decision**: Compute `amplitude[3]` for all event-observer families (threshold, Schmitt, neighborhood, and local-extremum), measured in `fVarIx` with the same semantics.

**Implementation Approach**: Phase 3 rolls out amplitude tracking to all semantic event families using the same extrema bookkeeping pattern.

### 5. Schmitt-Specific Readouts (Up/Down Durations, Duty Cycle) Are Family-Local

**Strategic Finding**: These statistics are specific to the Schmitt state machine and not generalizable to other trigger types.

**Strategic Decision**: Keep these as Schmitt-family-specific. Semantic Schmitt families may compute these in future work for feature parity with legacy `normalized_schmitt_trigger`, but they are not required for all event detectors.

**Implementation Approach**: Currently only legacy `normalized_schmitt_trigger` computes these. Defer semantic Schmitt family enhancement until phase 3-4 when period/amplitude decisions stabilize.

### 6. `min_amp` Semantics Are Defensible and Context-Dependent

**Strategic Finding**:
- 1-pass families: `min_amp` acts as a live-trajectory range gate
- 2-pass families: `min_amp` is compared against warmup-derived amplitude

This is not a bug; it's a consequence of the detection strategy.

**Strategic Decision**: Document clearly in user docs. This is defensible; users selecting 2-pass observers should expect 2-pass behavior.

### 7. Standardize on 2- or 3-Step Buffers with Best Helpers

**Strategic Decision**: Use consistent time, state, AND slope buffers across all families (not just state/time). Employ the best available:
- Compensated-time tracking for `elapsedTotal` and time-integrated means
- Best interpolation routines for event-time refinement (already using bounded three-sample quadratic for extrema)
- Shared helpers: `advanceAcceptedStepHistory2`, `advanceAcceptedStepHistory3`, `runningMeanTime`, `localMaximumFromThreeSamples`, etc.

**Benefit**: Consistent event-time accuracy across families and cleaner kernel code.

### 8. Active Dip Is Too Niche for Standardization

**Strategic Finding**: Only `normalized_schmitt_trigger` computes `activeDip`; it is specific to Schmitt down-state dynamics.

**Strategic Decision**: Keep as `normalized_schmitt_trigger`-specific legacy feature. Do not generalize.

### 9. Event State Value Selection Is a Separate Future Concern

**Strategic Finding**: All current observers store only event times, not event state or auxiliary values. Future goal is to return configurable event state subsets (like `solve_ivp`'s `dense_output`).

**Strategic Decision**: This is a separate architectural concern requiring a different seam (per-event state buffer vs. aggregate statistics). Track in Phase 6; do not implement yet.

### 10. Local-Maximum-First Unification Is A Viable Simplification Path

**Strategic Finding**: If oscillation statistics are the primary target, `local_maximum` can serve as a canonical base flow with optional storage of max/min event streams.

**Strategic Decision**: Use this as the default simplification direction for near-term observer-readout work. Keep `local_max` compatibility behavior stable, but avoid new local-extremum-specific surface growth until we decide whether it is retained as first-class or aliased/deprecated in favor of the canonical `local_maximum` contract.

### 11. Helper Evaluation and Interpolation Cleanup Should Be A Separate Backlog Stream

**Strategic Decision**: A rigorous helper review (compensated-time accumulation, compensated-mean primitives, interpolation variants such as robust 3-point Hermite/quadratic choices) is valuable but does not block current Phase 2-3 readout rollout.

**Implementation Approach**: Track this as helper-hardening backlog work, with dedicated numerical tests and benchmarking.

---

## Implementation Plan

### Phase 1: Audit & Documentation (✅ Complete)

1. ✅ Inventory current readouts across all families
2. ✅ Identify strategic decisions based on trigger geometry and ODE theory
3. ✅ Draft architectural recommendations with clear rationale
4. ✅ Clarify that semantic families are intentionally lean shells
5. ✅ Document context-dependent decisions (maxima counting, period naming, amplitude scope)

### Phase 2: Universal Readout Addition (Trajectory Stats + Period + Maxima Count)

**Scope**: Add readouts uniformly to all semantic families based on strategic decisions.

**Families to update**:
- `threshold_crossing`, `normalized_threshold_crossing`: Add trajectory stats, period, maxima count
- `schmitt_trigger`, `normalized_schmitt_trigger`: Add trajectory stats, period, maxima count
- `local_max`: Add trajectory stats, period (maxima count intentionally omitted)
- `normalized_neighborhood_return`: Add trajectory stats, period (maxima count validation needed separately)

**Readouts to add** (to all semantic families except as noted):
1. **Trajectory extrema** (universal): `max {var}`, `min {var}` per state variable
2. **Trajectory derivative extrema** (universal): `max d{var}/dt`, `min d{var}/dt` per state variable
3. **Trajectory mean** (universal): `mean {var}` per state variable (time-integrated with compensated time)
4. **Auxiliary trajectory statistics** (universal): `max {aux}`, `min {aux}`, `mean {aux}` per auxiliary variable
5. **Period** (universal): `max period`, `min period`, `mean period` (inter-event intervals, using elapsed time)
6. **Maxima count** (all except `local_max`): `max peak count`, `min peak count`, `mean peak count`

**Implementation approach**:
- Use best available compensated-time helpers for elapsed-time tracking
- Use bounded three-sample quadratic refinement for event-time accuracy
- Use `advanceAcceptedStepHistory2/3` and variants for buffer management
- Use `runningMean()` and `runningMeanTime()` for running aggregation

**Code changes**:
- `clode/kernels/observers/observer_threshold_crossing.clh`, etc.: Add trajectory-stats, period, maxima-count fields and updates
- `clode/observers/_definitions.py`: Update feature-name lists to include new readouts
- Tests: Validate new readouts per family

**Timeline**: 1-2 PRs (grouped by family type)

### Phase 3: Universal Amplitude Rollout Across Event Observers

**Scope**: Add amplitude statistics consistently to semantic event observers while converging oscillation semantics on the `local_maximum` contract.

**Readouts to add**:
1. **Amplitude**: `max amplitude`, `min amplitude`, `mean amplitude` for `threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, `local_max`, and `normalized_neighborhood_return`
2. **Retain existing period/maxima decisions** from Phase 2 (including no maxima count for `local_max`)
3. **Compatibility rule**: do not add new local-extremum-only readout concepts while this phase is in flight

**Implementation approach**:
- Track last-seen local maxima and minima in `fVarIx` for each family
- Compute amplitude as `last_local_max - last_local_min` on valid updates
- Use existing three-sample refinement and compensated-time history buffers
- Reuse one oscillation-statistics helper pattern across families to avoid duplicating policy logic

**Code changes**:
- `clode/kernels/observers/observer_*.clh` (semantic event families): Add amplitude tracking fields and updates
- `clode/observers/_definitions.py`: Update feature-name lists for amplitude outputs
- Tests: Validate amplitude consistency against `local_maximum` on representative oscillators

**Validation**: Separate exact-solution or limit-cycle comparison with `local_maximum` as canonical oscillation reference.

**Timeline**: 1 PR

### Phase 4: Family-Alignment Decision (Schmitt Extras + Simplification Options)

**Scope**: Decide whether remaining gaps (Schmitt-specific readouts, local-extremum compatibility cleanup) warrant further expansion.

**Questions to resolve**:
- Should semantic Schmitt families compute up/down durations and duty cycle for feature parity with `normalized_schmitt_trigger`?
- Should `local_max` remain first-class, alias to `local_maximum` semantics, or move to compatibility/deprecation path?

**Decision output**: Update design notes; may proceed to partial Phase 5 or defer extremum/neighborhood enhancement to future work.

**Timeline**: Planning only; no implementation

### Phase 5: Selective Readout Surface (Conditional on Phase 4)

**Scope**: Decide whether to introduce family-level readout-selection surfaces at all (none are currently public).

**Conditional implementation** (if strategic alignment is achieved in Phase 4):
- Consolidate helper functions in `clODE_observer_helpers.clh` per audit recommendations
- Implement readout-selection types for families with optional readouts only where evidence supports the added surface
- Extend selective surfaces to threshold, Schmitt, and extremum families

**Timeline**: 2-3 PRs (if Phase 4 approves) or documentation only (if family-specific is final decision)

---

## Implementation Checklist and Next Actions

**For Phase 2 (Next PR)**:
- [x] Add trajectory-stats fields to `observer_threshold_crossing.clh`, `observer_normalized_threshold_crossing.clh`
- [x] Add period computation using elapsed time and running means
- [x] Add maxima-count computation in these families
- [x] Update `clode/observers/_definitions.py` feature-name resolution
- [x] Update tests and documentation
- [x] Validate no behavior change to event detection; only new output schema

**For Phase 3 (Separate PR)**:
- [ ] Add amplitude tracking to all semantic event observers
- [ ] Add opposite-polarity extrema bookkeeping for robust amplitude updates
- [ ] Update feature-name resolution for amplitude outputs
- [ ] Validate output consistency with `local_maximum` on test models
- [ ] Decide and document `local_max` compatibility path (`retain` vs `alias` vs `deprecate`)

**Ongoing**:
- [ ] Keep readout-selection surface decisions evidence-driven; do not assume a selector type is required for every family
- [ ] Document context-dependent meanings of "period" and "amplitude" per family
- [ ] Reference audit findings in docs to justify design decisions
- [ ] Add helper-hardening backlog item for compensated-time/mean and interpolation routine review
