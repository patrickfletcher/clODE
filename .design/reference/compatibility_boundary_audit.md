# Compatibility Boundary Audit

Purpose: record which observer-related public surfaces are canonical, which remain compatibility-only, and which names should stay internal-only.
Read when: changing observer naming, config surfaces, root exports, or public docs for canonical observer families.
Update when: a compatibility alias is added, removed, deprecated, promoted to canonical status, or explicitly limited to compatibility-only documentation.

## Current boundary

### Canonical semantic surfaces

- `Observer.threshold_crossing`, `Observer.normalized_threshold_crossing`, `Observer.schmitt_trigger`, and `Observer.normalized_schmitt_trigger` are the preferred public observer names for the current threshold catalog.
- `ThresholdCrossingConfig` is the preferred semantic config surface for the two one-boundary threshold observers.
- `SchmittTriggerConfig` is the preferred semantic config surface for the two canonical Schmitt families.
- `Observer.local_max` is the preferred public observer name for extremum workflows, and `LocalMaximumConfig` is its preferred semantic config surface.
- `Observer.normalized_neighborhood_return` is the preferred public observer name for the canonical two-pass normalized neighborhood-return family, and `NeighborhoodReturnConfig` is its preferred semantic config surface.
- The selected observer, not the config class alone, determines whether a shared config is interpreted in absolute units or as warmup-derived fractions.

### Supported compatibility surfaces

- `Observer.normalized_schmitt_trigger` keeps derivative-gate controls on the compatibility surface; those controls are no longer part of `SchmittTriggerConfig`.
- Observer names in the `Observer` enum are canonical; compatibility is provided via parameter/config surfaces and constructor keyword adapters.
- `ObserverParams` remains the broad compatibility bundle for observer settings, including canonical event-trigger controls such as `max_event_count` and `min_amp`, plus compatibility-only fields such as `eps_dx` and `min_imi`.
- `FeatureSimulator` constructor `observer_*` keyword arguments and `set_observer_parameters(...)` remain supported compatibility paths.
- Flat root barrels and re-exports such as `clode.__init__`, `clode.features`, and the older top-level compatibility modules remain compatibility import surfaces unless a later audit explicitly promotes or deprecates them.

### Internal-only surfaces

- Raw observer implementation names such as `thresh1`, `thresh2`, `thresh3`, and `schmitt` are internal build identifiers, not public API.
- Kernel/header filenames such as `observer_threshold_crossing.clh`, `observer_normalized_threshold_crossing.clh`, `observer_schmitt_trigger.clh`, `observer_normalized_schmitt_trigger.clh`, and `observer_normalized_neighborhood_return.clh` are implementation details.
- Normalized runtime-setting structs such as `ObserverRuntimeSettings` and `EventOutputSettings` are internal semantic owners, not public compatibility promises.

## Current decisions

- Removed historical observer labels are not part of the current public observer surface.
- Observer enum names are canonical for public use.
- User docs should prefer semantic names and shared config classes, and should mention compatibility aliases only when that helps users translate older code.
- Public docs should describe `ObserverParams` and `observer_*` keyword arguments as compatibility surfaces rather than as the preferred observer UX.
- Future observer-family work should treat the compatibility bundle as an adapter boundary, not as the vocabulary source for new semantic config objects.
- `SchmittTriggerConfig` round-trips both semantic Schmitt variants, including `Observer.normalized_schmitt_trigger`, through `get_observer_configuration()` and `set_observer_configuration(...)`.
- Schmitt state-machine duration outputs (`up duration`, `down duration`, `duty`) are considered core Schmitt measurements for canonical observer behavior.
- Threshold, Schmitt, and neighborhood semantic configs now expose `feature_var` where their kernels split trigger geometry from extrema/amplitude semantics.
- `max_event_count` is a canonical event-trigger-family limiter and early-stop control, even though it still passes through `ObserverParams` on compatibility paths.
- `min_amp` should be documented as a user-specified gate on variation in `event_var`; it is now canonical across the current event-triggering semantic families, with one-pass and warmup-seeded families differing only in how they accumulate that variation.
- `eps_dx` has no current semantic-family use case and remains compatibility-only pending a dedicated audit.
- The readout/kernel audit did not justify a new cross-family "oscillation bundle" config object. Keep family-specific semantic config objects, and treat future sharing as build-specialized readout bundle declarations or selectors instead.

## Questions for follow-through

- Should `ObserverParams` and the `observer_*` compatibility paths eventually receive either a formal deprecation posture or a tighter long-term support statement once family-specific config surfaces settle further?
- When the package eventually audits flat compatibility barrels, which ones should remain public convenience surfaces and which should move toward retirement?
- How should build-specialized readout bundle declarations or selectors surface without re-expanding `ObserverParams` or the `observer_*` compatibility paths?

## Touchpoints

- `clode/observers/types.py`
- `clode/observers/_definitions.py`
- `clode/simulation/features.py`
- `clode/__init__.py`
- `clode/observers/__init__.py`
- `docs/feature_extraction.md`
- `docs/observers.md`
- `.design/next_pr.md`
- `.design/package_state.md`
