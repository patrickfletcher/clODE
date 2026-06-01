# Observer Readout Architecture Phases

Purpose: track phased observer-readout decisions and what remains open after landed semantic observer rollouts.
Read when: evaluating follow-on observer readout parity, family-specific extras, or potential future readout-selection surfaces.
Update when: a readout phase lands, is deferred, or changes scope materially.

This document lays out the planned phases for expanding observer readout-selection surfaces and resolving family-alignment questions. Each phase has explicit dependencies, success criteria, and blocking points.

**Overall Goal**: Enable users to choose which readouts they want from each observer family (narrow output schema) without having to understand kernel-private storage details or architectural seams.

**Current State** (after latest semantic rollout):
- Phase 1 ✅: Audit complete; readout categorization and architecture recommendations documented
- Phase 2 ✅: Core trajectory/period/maxima rollout is implemented and regression-validated in the active PR slice
- Phase 3 ✅: Amplitude rollout is implemented across semantic event families and validated on focused observer-contract regressions
- Phase 4: Decision-focused planning remains open (see `.design/next_pr.md`)
- Phase 5: Conditional future work; no public readout-selection dataclasses are currently live

---

## Phase 1: Audit and Architecture Alignment (✅ Complete)

**Scope**: Comprehensive audit of readouts across all observer families.

**Deliverables**:
- ✅ `.design/reference/observer_readout_audit.md` — Full categorization, inventory, findings, and recommendations
- ✅ `.design/package_state.md` updated — Observer-readout architecture summary
- ✅ `.design/next_pr.md` updated — PR scope and acceptance criteria
- ✅ `docs/observers.md` updated — Readout inventory table and `min_amp` semantics clarified
- ✅ `docs/feature_extraction.md` updated — semantic observer workflows and readout behavior clarified

**Success Criteria**:
- ✅ Clear categorization of readouts (universal, family-specific, cross-family-not-yet-aligned)
- ✅ Phase recommendations documented with blocking dependencies
- ✅ User docs clarified without historical narration
- ✅ Design notes agree on next priority

**Blocking Dependencies**: None (phase 1 is foundational)

---

## Phase 2: Trajectory Statistics + Period + Maxima Count Alignment

**Goal**: Add readouts uniformly to all semantic families based on strategic decisions; establish baseline feature parity across trigger types.

**Strategic Context**: Semantic families are lean shells; Phase 2 adds foundational readouts to all of them (trajectory stats, period, maxima count). This is **not** a readout-selection feature; it fills gaps and establishes consistent coverage across families.

**Rationale**:
- **Trajectory statistics** (extrema, means) are independent of event-detection mechanism and are universally useful
- **Period** ("inter-event interval") is meaningful for all event detectors with clear, consistent computation
- **Maxima counting** is used for non-extremum semantic triggers and intentionally omitted for `local_max` to reduce overlap with the canonical `local_maximum` oscillation contract

**Scope**:
- Add `xTrajectoryMax[N_VAR]`, `xTrajectoryMin[N_VAR]`, `xTrajectoryMean[N_VAR]` to all semantic families
- Add `dxTrajectoryMax[N_VAR]`, `dxTrajectoryMin[N_VAR]` to all semantic families
- Add `auxTrajectoryMax[N_AUX]`, `auxTrajectoryMin[N_AUX]`, `auxTrajectoryMean[N_AUX]` to all semantic families
- Add `period[3]` (max/min/mean inter-event interval, using elapsed time) to all semantic families
- Add `nMaxima[3]` (max/min/mean local maxima count) to all semantic families **except** `local_max`

**Families Affected**:
- `threshold_crossing`, `normalized_threshold_crossing`
- `schmitt_trigger`, `normalized_schmitt_trigger`
- `local_max` (no maxima count)
- `normalized_neighborhood_return`

**Implementation Notes**:
- Use **compensated elapsed time** (not wallclock `t`) for period and time-integrated means
- Use bounded three-sample quadratic refinement (already available) for event-time accuracy
- Leverage `runningMeanTime()` and `advanceAcceptedStepHistory3()` helpers
- Maxima detection: track sign changes in derivative (already in use)

**State Overhead**:
- Per-family addition: `5*N_VAR + 3*N_AUX + 3*3` floats (period, nMaxima are each 3-element arrays)
- Total: modest (~50-100 floats for typical N_VAR=10, N_AUX=5)

**Code Changes**:
- `clode/kernels/observers/observer_threshold_crossing.clh`, `observer_normalized_threshold_crossing.clh`: Add fields and tracking
- `clode/kernels/observers/observer_schmitt_trigger.clh`, `observer_normalized_schmitt_trigger.clh`: Add fields and tracking
- `clode/kernels/observers/observer_local_maximum.clh`: Add trajectory stats, period (NOT maxima count)
- `clode/kernels/observers/observer_normalized_neighborhood_return.clh`: Add trajectory stats, period, maxima count
- `clode/observers/_definitions.py`: Update feature-name lists to include new readouts
- Tests: Validate new readouts per family; no behavior change to event detection

**Acceptance Criteria**:
- All semantic families emit trajectory statistics with same names and semantics as legacy families
- Period is computed and named consistently across all families
- Maxima counting is present for threshold, Schmitt, neighborhood families; absent for local-extremum
- Regression tests pass: existing event-detection behavior unchanged
- New tests validate trajectory stats are populated correctly
- Documentation updated to clarify readout availability and meaning per family

**Blocking Dependencies**: Phase 1 ✅

**Timeline Estimate**: 1-2 PRs (grouped by family type or combined)

---

## Phase 3: Universal Amplitude Rollout With Local-Maximum-First Simplification (✅ Implemented)

**Goal**: Add amplitude statistics to semantic event observers while converging on one canonical oscillation contract (`local_maximum` semantics).

**Strategic Context**: After Phase 2, all semantic families have trajectory stats and period, and most have maxima counting. Phase 3 adds amplitude tracking as a universal oscillation statistic for event observers, computed in `fVarIx` via local extrema bookkeeping independent of trigger geometry. To reduce drift and bloat, no new local-extremum-only readout behavior should be introduced while this phase lands.

**Delivered scope**:
- Add amplitude computation to `threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, `local_max`, and `normalized_neighborhood_return`
- Add `amplitude[3]` (max/min/mean) to semantic family storage and output schema
- Keep `period[3]` semantics from Phase 2 unchanged
- **Intentional omission remains**: `nMaxima[3]` is NOT computed for `local_max`
- Compatibility decision for `local_max` is deferred to Phase 4

**Code changes**:
- `clode/kernels/observers/observer_*.clh` (semantic event families): Add opposite-polarity extrema tracking and amplitude computation
- `clode/observers/_definitions.py`: Update feature-name lists to include amplitude across semantic event families
- Tests: Validate amplitude matches `local_maximum` on representative oscillation models

**Validation status**:
- Focused contract regressions pass after rollout (`test/test_features.py` and `test/test_simulation_contracts.py`)
- Direct `local_maximum` parity benchmarking remains follow-on evidence for Phase 4 decision support

**State Overhead**:
- Additional per family: extrema-tracking fields plus `amplitude[3]`

**Acceptance criteria status**:
- ✅ All semantic event families emit `max/min/mean amplitude`
- ✅ Maxima counting remains absent for `local_max` (by simplification policy)
- ✅ Focused regressions pass with unchanged event-detection contracts
- ⏳ Remaining: document final compatibility direction for `local_max`

**Blocking Dependencies**: Phase 2 ✅ (for trajectory stats and period; Phase 3 built on them)

**Timeline Estimate**: 1 PR

---

## Phase 4: Family-Alignment and Compatibility Cleanup Decision (Point of Commitment, Active)

**Goal**: Decide how to handle edge cases and remaining readout possibilities after phases 2-3 establish baseline coverage.

**Strategic Context**: After phases 2-3, we have:
- Universal trajectory stats, period, and maxima counting (except local-extremum)
- Enhanced `local_max` with amplitude (no maxima counting by design)
- Clear architectural seams for core readouts

Phase 4 addresses remaining questions about edge cases and optional enhancements:

**Questions to Decide**:
1. Should semantic Schmitt families compute up/down durations and duty cycle (for feature parity with legacy `normalized_schmitt_trigger`)?
   - Current: only legacy `normalized_schmitt_trigger` computes these
   - Option A: Add to semantic Schmitt in Phase 5
   - Option B: Keep as legacy-only; document as such
2. Should `local_max` remain first-class, alias to canonical `local_maximum` oscillation semantics, or move to a compatibility/deprecation path?
   - Current: both surfaces coexist and can drift
   - Question: which option preserves clarity and minimizes maintenance cost?
3. Are there other edge cases or family-specific readouts not yet covered?

**Scope**:
- Review phases 2-3 outcomes: what works, what diverges
- Assess user need for Schmitt-specific readouts in semantic families
- Validate neighborhood-family readout semantics
- Document decisions and rationale

**Code Changes**: None (decision-only phase; implementation follows decisions in Phase 5)

**Deliverables**:
- Updated `.design/reference/observer_readout_audit.md` with decision and rationale
- Updated `.design/next_pr.md` with follow-on plan
- If decisions warrant implementation: add recommendations to Phase 5 scope

**Acceptance Criteria**:
- Clear decision on Schmitt-specific readout expansion: Phase 5 scope or future work
- Clear decision on neighborhood-family readout validation: part of Phase 2 validation or separate
- Rationale documented in design notes
- No code changes; decisions guide Phase 5

**Blocking Dependencies**: Phases 2 and 3 ✅

**Timeline Estimate**: Planning and documentation only (~1 day)

---

## Phase 5: Selective-Readout Surface and Helper Consolidation (Conditional)

**Goal**: define whether family-level readout-selection surfaces are worth adding, and consolidate helper functions where duplication is measurable.

**Strategic Context**: After phases 1-4, architectural seams are clear and core readouts are implemented. Phase 5 addresses the user-facing question: "Can I choose which readouts I want from each family?"

**Conditional Implementation**:

**If Phase 4 approves Schmitt-specific enhancement**:
- Implement up/down duration and duty cycle computation in semantic Schmitt families
- Extend period/amplitude/maxima-count to any families that lacked them in phases 2-3

**Core Work** (regardless of Phase 4 decision):
1. **Consolidate helpers** in new or updated `clODE_observer_helpers.clh`:
   - Formalize `updateTrajectoryStats()`, `updatePeriodStats()`, `updateAmplitudeStats()`, `updateMaximaCountStats()`
   - Leverage compensated-time and best interpolation routines
   - Reduce code duplication across `observer_threshold_crossing.clh`, `observer_schmitt_trigger.clh`, `observer_local_maximum.clh`, `observer_normalized_neighborhood_return.clh`

2. **Optionally implement readout-selection types** for families with optional readouts:
   - no readout-selection dataclass is currently part of the public observer surface
   - any future selector types should be introduced only after Phase 4 decisions are explicit

3. **Update feature-schema resolution** to handle selective readouts at build time

**Code Changes** (if proceeding):
- `clode/kernels/clODE_observer_helpers.clh`: Formalized and consolidated helpers
- `clode/observers/types.py`: New readout-selection dataclasses
- `clode/observers/_definitions.py`: Resolution logic for new selectors
- `clode/_opencl/executors.py`: Build and metadata routing for new selectors
- Kernel definitions: Conditional compilation guards for optional readouts
- Tests: Coverage for new selector combinations per family

**Acceptance Criteria** (if proceeding):
- All applicable families expose readout-selection surfaces (or explicit rationale for why they don't)
- Naming is short and consistent across families; meaning disambiguated by documentation
- Regression tests pass: existing behavior unchanged
- New tests cover selector combinations per family
- Documentation updated with new selector types and examples
- Performance: no overhead when readouts are not selectively disabled
- Helper consolidation reduces kernel code duplication

**If Phase 4 Decision Is "Family-Specific" or "Defer"**:
- Document explicit rationale for why selective readouts remain family-scoped
- Clarify which families support selection and which do not
- Consider Phase 6 (event-state values) or alternative priorities

**Blocking Dependencies**: Phase 4 decision ✅

**Timeline Estimate**: 2-3 PRs (if helper consolidation + multiple selector types) or documentation only (if family-specific/defer decision)

**Acceptance Criteria** (if universal decision):
- All applicable families expose readout-selection surfaces
- Naming is short and consistent across families (context disambiguates)
- Regression tests pass: existing behavior unchanged
- New tests cover selector combinations per family
- Documentation updated with new selector types and examples
- Performance: no overhead when readouts are not selectively disabled

**Acceptance Criteria** (if family-specific decision):
- Each family's readout complement clearly documented
- Rationale for family-specific readouts stated in design notes
- Next alternative priority identified and added to `.design/ideas.md`

**Blocking Dependencies**: Phase 4 decision ✅

**Blocked By**: Phase 4 decision

**Timeline Estimate**: 2-3 PRs (if universal decision); documentation only (if family-specific decision)

---

## Phase 6: Event State Values (Future, Lower Priority)

**Goal**: Design and implement configurable per-event state capture (e.g., all state variables at each event, not just times).

**Rationale**: `solve_ivp`-style behavior is to return not just event times but also state values at events. This is a separate architectural concern from aggregate statistics and should not block phases 2-5.

**Scope**:
- Design configurable per-event state capture (which variables, which events)
- Implement storage and retrieval mechanisms
- Add to event-detector families

**Acceptance Criteria**:
- Users can request full or partial state values at events
- Output schema reflects user selection
- No overhead for users not requesting event state values

**Blocking Dependencies**: Phases 2-4 conceptual clarity (no hard blocker)

**Timeline Estimate**: TBD (lower priority; may be deferred indefinitely if not demanded)

---

## Dependency Graph

```
Phase 1 ✅ (Audit)
  ↓
  ├→ Phase 2 (Trajectory Stats)
  │    ↓
  │    └→ Phase 4 (Family Alignment Decision)
  │         ↓
  │         └→ Phase 5 (Selective Readout UX)
  │
  └→ Phase 3 (Semantic Extremum)
       ↓
       └→ Phase 4 (Family Alignment Decision)
            ↓
            └→ Phase 5 (Selective Readout UX)

Phase 6 (Event State Values) — independent, lower priority
```

**Critical Path**: Phase 1 → Phase 2 + Phase 3 (parallel) → Phase 4 → Phase 5 (optional)

---

## Decision Points and Rollback Paths

### Phase 4 Decision Point

After phases 2-3, the following **already decided** by strategic audit (not open decisions):
- ✅ Trajectory statistics are universal (all families)
- ✅ Period ("inter-event interval") is universal (all families, context-dependent meaning)
- ✅ Maxima counting is implemented on non-extremum semantic families; `local_max` omits it by simplification policy
- ✅ Amplitude is universal across event observers (computed in `fVarIx` via extrema tracking)
- ✅ `local_maximum` semantics are the canonical oscillation reference model

**Open Decision**: Additional optional enhancements beyond core Phase 2-3 implementation:

**If Phase 4 approves additional Schmitt enhancement**:
- Compute up/down durations and duty cycle in semantic Schmitt families (for feature parity with legacy `normalized_schmitt_trigger`)
- Proceed to Phase 5 with broader selective-readout expansion

**If Phase 4 decides defer Schmitt enhancement**:
- Keep as legacy-only feature
- Proceed to Phase 5 with core helper consolidation but narrower selectivity

**If any edge cases or validation issues surface**:
- Document explicitly; may delay Phase 5 or narrow its scope

**Outcome**: Updated design notes with approval for Phase 5 scope (either full or partial) or rationale for deferral

---

## Success Metrics (End-to-End)

**After Phases 1-3** (core readouts implemented):
- All semantic event-detector families emit trajectory statistics (max/min/mean per variable)
- All families compute period as "inter-event interval" with consistent naming
- Non-extremum semantic families compute maxima count per period; `local_max` omits it by simplification policy
- All semantic event families compute amplitude (opposite-polarity extrema tracking)
- `min_amp` semantics documented clearly per family type
- Regression tests pass: existing event-detection behavior unchanged
- No performance regression; modest state overhead (~50-100 floats per instance)

**After Phases 1-4** (decision point):
- Phase 4 decision documented and rationale clear
- Any Schmitt-specific enhancement decision made
- Phase 5 scope finalized (full helper consolidation + selectors, or partial, or deferred)

**After Phases 1-5** (selective-readout surfaces, if Phase 4 approves):
- Users can inspect available readouts per family via `get_feature_names()`
- Users can choose readout subsets where applicable (e.g., period only, amplitude only, etc.)
- Helper consolidation in `clODE_observer_helpers.clh` reduces kernel code duplication
- Performance: no overhead when readouts are not selectively disabled
- Short, context-dependent naming used consistently across families

**Across all phases**:
- No architectural debt introduced
- Readout selection remains narrow and well-justified
- Documentation is current and avoids historical narration
- Semantic families achieve feature parity with legacy where appropriate

---

## Notes and Strategic Context

### Already Decided (From Phase 1 Audit):
- Semantic families are **intentionally lean shells** designed to have readouts added consistently in phases 2-3
- "Period" = "inter-event interval" universally; meaning is context-dependent by trigger geometry
- Maxima counting is kept on non-extremum semantic families; `local_max` omits it to avoid duplicating the canonical `local_maximum` oscillation contract
- Amplitude should be available for **all event observers** with consistent semantics
- `min_amp` semantics are defensible but context-dependent (1-pass vs. 2-pass); documentation is sufficient
- Use compensated elapsed time, bounded three-sample refinement, and best available helpers for consistency

### Conditional / Open (For Phase 4 Decision):
- Whether semantic Schmitt families should compute up/down durations and duty cycle (for feature parity with legacy `normalized_schmitt_trigger`)
- Whether `local_max` should be retained, aliased, or deprecated as a compatibility surface
- Whether to consolidate helpers in Phase 5 or defer
- Whether selective-readout surfaces should be introduced at all (none are currently live)

### Will Not Change:
- Event-state-value selection is a separate architectural concern (Phase 6, much later)
- Any future readout-selection surface should remain optional and evidence-driven rather than a blanket requirement across families
- Core Phase 2-3 scope is fixed: trajectory stats, period, maxima count (except local-extremum), amplitude (all event families)
