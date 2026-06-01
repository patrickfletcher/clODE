# Observer Readout Audit

Purpose: capture the current observer readout contract, resolved policy decisions, and evidence-backed guardrails for follow-on observer work.
Read when: deciding observer readout behavior, naming, or schema boundaries across threshold, Schmitt, local-extremum, and neighborhood-return families.
Update when: a readout schema changes, a family gains/removes readouts, or an open decision in this note is resolved.

## Fast path

- Use this note for current readout facts and open decisions.
- Use `.design/next_pr.md` for active implementation scope and acceptance criteria.
- Use `.design/reference/observer_concept_audit.md` for deeper architecture rationale.
- Use `.design/reference/compatibility_boundary_audit.md` for canonical vs compatibility surface policy.
- Historical phased rollout material is archived in `.design/archived/observer_readout_cleanup_2026_05_31/observer_readout_audit.md`.

## Current readout contract

### Summary family

- `Observer.summary` emits selected summary reductions (state, auxiliary, and slope groups) from `SummaryObserverSelection`.
- No event streams are emitted for this family.

### Event-observer families

All semantic event observers now emit event timestamps plus observer-local summary/readout groups from `clode/observers/_definitions.py`.

| Family | Event streams | Oscillation-oriented readouts | Notes |
| --- | --- | --- | --- |
| `Observer.threshold_crossing` | one event-time stream + event count | period, maxima count, amplitude | one-pass absolute threshold semantics |
| `Observer.normalized_threshold_crossing` | one event-time stream + event count | period, maxima count, amplitude | two-pass warmup-derived threshold semantics |
| `Observer.schmitt_trigger` | up/down transition streams + event count | period, maxima count, up/down duration, duty, active dip, amplitude | one-pass absolute Schmitt semantics |
| `Observer.normalized_schmitt_trigger` | up/down transition streams + event count | period, maxima count, up/down duration, duty, active dip, amplitude | two-pass warmup-derived Schmitt semantics |
| `Observer.local_max` | max/min event time + value streams + event count | IMI-period statistics, amplitude | intentionally omits maxima count |
| `Observer.normalized_neighborhood_return` | one event-time stream + event count | period, maxima count, amplitude | includes neighborhood range/anchor-oriented outputs |

### Cross-family naming contract

- For threshold/Schmitt/neighborhood families, `period` names mean inter-event interval; interpretation is family-specific by trigger geometry.
- For `local_max`, inter-event interval is surfaced as IMI (`max IMI`, `min IMI`, `mean IMI`) in the current schema.
- `amplitude` is tracked from extrema bookkeeping in the family-specific feature channel. Threshold, Schmitt, and neighborhood-return use `f_var_ix`; `local_max` uses its one resolved extrema channel.
- `n maxima` is present for non-extremum event families and intentionally absent for `local_max`.
- Schmitt-only state-machine extras now include `up duration`, `down duration`, `duty`, and the family-local `active dip` readout.

### Variable-routing implications

- Threshold event detection is keyed to `e_var_ix`, while threshold-family extrema and amplitude tracking now use `f_var_ix`.
- Schmitt event detection is keyed to `e_var_ix`, while maxima/amplitude/active-dip tracking uses `f_var_ix`.
- `local_max` detection is keyed to `f_var_ix`.
- Neighborhood-return uses `e_var_ix` for anchor/threshold behavior and `f_var_ix` for maxima/amplitude tracking.
- See `.design/reference/observer_event_feature_variable_contract.md` for the full routing contract.

## Resolved policy decisions

- Schmitt-specific duration outputs (`up duration`, `down duration`, `duty`) are core Schmitt measurements and are now exposed on canonical Schmitt readout surfaces.
- Threshold, Schmitt, and neighborhood semantic configs now expose `feature_var` wherever their live kernels split trigger geometry from extrema/amplitude behavior; `local_max` remains the one-channel exception.

## Open policy decisions

- Whether recurring event-trigger controls plus oscillation-oriented readouts should be factored through a shared config seam or remain family-local: active target in `.design/next_pr.md`.
- Whether a family-level readout-selection surface should be introduced: an open packaging question; any such surface should not depend on historical family labeling.

## Cross-family control semantics

- `max_event_count` is the current cross-family event-loop limiter across the event-triggering observer families; users can treat it as an early-stop control when they only need a bounded number of events.
- `min_amp` is a user-specified gate on variation in `event_var` that suppresses tiny or noisy oscillations when they are not of interest. It is now part of every current event-triggering semantic config surface, with one-pass and warmup-seeded families differing only in how that variation is accumulated.
- `eps_dx` and `min_imi` remain compatibility-only controls on `ObserverParams`; they do not currently define the preferred semantic observer UX.

## Guardrails

- Prefer one canonical semantic meaning per public readout name; avoid reintroducing legacy aliases in new docs.
- Treat `ObserverParams` and constructor `observer_*` settings as compatibility adapters, not as the source of new semantic readout vocabulary.
- Keep rollout sequencing and implementation history out of this live note once work lands; archive historical planning detail instead.
- When a readout contract changes, update this note and the relevant user docs in the same pass.

## Evidence anchors

- `clode/observers/_definitions.py`
- `clode/observers/types.py`
- `clode/simulation/features.py`
- `clode/kernels/observers/observer_threshold_crossing.clh`
- `clode/kernels/observers/observer_schmitt_trigger.clh`
- `clode/kernels/observers/observer_local_maximum.clh`
- `clode/kernels/observers/observer_normalized_neighborhood_return.clh`
- `test/test_features.py`
- `test/test_simulation_contracts.py`
