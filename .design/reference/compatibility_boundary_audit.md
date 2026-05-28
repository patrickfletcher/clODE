# Compatibility Boundary Audit

Purpose: record which observer-related public surfaces are canonical, which remain compatibility-only, and which names should stay internal-only.
Read when: changing observer naming, config surfaces, root exports, or public docs for semantic observer families and retained legacy observers.
Update when: a compatibility alias is added, removed, deprecated, promoted to canonical status, or explicitly limited to compatibility-only documentation.

## Current boundary

### Canonical semantic surfaces

- `Observer.threshold_crossing`, `Observer.normalized_threshold_crossing`, `Observer.schmitt_trigger`, and `Observer.normalized_schmitt_trigger` are the preferred public observer names for the current threshold catalog.
- `ThresholdCrossingConfig` is the preferred semantic config surface for the two one-boundary threshold observers.
- `SchmittTriggerConfig` is the preferred semantic config surface for the two lean Schmitt families.
- `Observer.local_extremum` is the preferred public observer name for the lean polarity-selectable extremum family, and `LocalExtremumConfig` is its preferred semantic config surface.
- `Observer.neighborhood_return` is the preferred public observer name for the lean two-pass normalized neighborhood-return family, and `NeighborhoodReturnConfig` is its preferred semantic config surface.
- The selected observer, not the config class alone, determines whether a shared config is interpreted in absolute units or as warmup-derived fractions.

### Supported compatibility and legacy surfaces

- `Observer.threshold_2` remains a supported legacy-named public observer for the fully featured warmup-derived normalized Schmitt workflow.
- `Observer.threshold_2` keeps derivative-gate controls on the compatibility surface; those controls are no longer part of `SchmittTriggerConfig`.
- `Observer.local_max` remains a supported legacy public observer for the heavier maxima-oriented extremum workflow.
- `Observer.neighbourhood_2` remains a supported legacy public observer for the heavier normalized neighborhood-return workflow.
- `ObserverParams` remains the broad compatibility bundle for observer settings.
- `FeatureSimulator` constructor `observer_*` keyword arguments and `set_observer_parameters(...)` remain supported compatibility paths.
- Flat root barrels and re-exports such as `clode.__init__`, `clode.features`, and the older top-level compatibility modules remain compatibility import surfaces unless a later audit explicitly promotes or deprecates them.

### Internal-only surfaces

- Raw observer implementation names such as `thresh1`, `thresh2`, `thresh3`, and `schmitt` are internal build identifiers, not public API.
- Kernel/header filenames such as `observer_threshold_crossing.clh`, `observer_normalized_threshold_crossing.clh`, `observer_schmitt_trigger.clh`, `observer_normalized_schmitt_trigger.clh`, and `observer_threshold_2.clh` are implementation details.
- Normalized runtime-setting structs such as `ObserverRuntimeSettings` and `EventOutputSettings` are internal semantic owners, not public compatibility promises.

## Current decisions

- There is no public `Observer.threshold_1` or `Observer.threshold_3` alias.
- `Observer.threshold_2` is a retained public legacy observer, not the semantic owner of normalized Schmitt naming.
- `Observer.local_max` is a retained public legacy observer, not the semantic owner of local-extremum naming.
- `Observer.neighbourhood_2` is a retained public legacy observer, not the semantic owner of neighborhood-return naming.
- User docs should prefer semantic names and shared config classes, and should mention compatibility aliases only when that helps users translate older code.
- Public docs should describe `ObserverParams` and `observer_*` keyword arguments as compatibility surfaces rather than as the preferred observer UX.
- Future observer-family work should treat the compatibility bundle as an adapter boundary, not as the vocabulary source for new semantic config objects.
- `SchmittTriggerConfig` no longer round-trips `Observer.threshold_2`; `get_observer_configuration()` is intentionally `None` for that retained legacy observer.

## Questions for follow-through

- Should `threshold_2` remain indefinitely as a retained public legacy observer, or only until a broader public API audit names a deprecation posture?
- Should `local_max` and `neighbourhood_2` remain indefinitely as retained public legacy observers, or only until the oscillation-oriented readout story is settled enough to name a deprecation posture?
- Should `ObserverParams` and the `observer_*` keyword path stay lightly documented compatibility helpers, or should they get a tighter long-term support statement?
- When the package eventually audits flat compatibility barrels, which ones should remain public convenience surfaces and which should move toward retirement?
- Does the current threshold-family reuse justify a shared oscillation-oriented bundle seam, or is the cleaner next step still to keep family-local config objects and duplication where needed?

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
