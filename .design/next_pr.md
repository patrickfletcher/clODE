# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Observer readout bundle/selectability proof

## Scope

Define the first shared readout-bundle vocabulary for the current event observers, keep family-specific semantic config objects in place, and land a build-specialized selection proof on the threshold-family observers.

Current audit findings that shape this target:

- The semantic event observers now have standardized kernel numerics and helper usage; the remaining work is packaging/selectability, not baseline mean/interpolation cleanup.
- The family-defining readouts are already close to maximal: threshold and neighborhood-return cover the one-stream oscillation core, Schmitt owns the phase-state extras, and `local_max` owns the extremum-stream/IMI path.
- The natural shared bundles are `event geometry`, `oscillation core`, and `trajectory summary`; Schmitt phase-state outputs and neighborhood/local-max extras should stay family-local for the first proof.
- `min_amp` and `max_event_count` already have semantic homes on family-specific config objects, so the next reuse seam should be on readout declarations rather than on another cross-family config bundle.

The proof target is: document the bundle inventory, choose the build-specialized selectability mechanism, and adopt it first on `Observer.threshold_crossing` plus `Observer.normalized_threshold_crossing` without changing their default public schemas.

**Strictly out of scope:**

- More kernel time/mean/interpolation refactors
- Observer-state or register-pressure cleanup
- New observer families
- A new cross-family semantic config object for recurring event-trigger controls

## Why Now

The helper-hardening pass is now landed: semantic event observers share K=3 history updates, compensated trajectory/auxiliary means, and family-consistent interpolation behavior. The readout audit also shows that the remaining ambiguity is not "which extra outputs are missing?" but "how should the existing outputs be grouped and selectively exposed without bloating every family?"

## Acceptance Criteria

- [ ] `.design/reference/observer_readout_audit.md` records the shared bundle inventory and the rationale for keeping family-specific semantic configs.
- [ ] `.design/reference/compatibility_boundary_audit.md` records that a control-only oscillation-bundle config is not the preferred next step.
- [ ] A build-specialized readout-selection surface exists for the threshold family with at least `event geometry`, `oscillation core`, and `trajectory summary` bundles.
- [ ] `Observer.threshold_crossing` and `Observer.normalized_threshold_crossing` adopt the new bundle declarations without changing their default public feature schemas.
- [ ] The family-local status of Schmitt phase-state outputs and neighborhood/local-max extras is documented explicitly.
- [ ] `package_state.md`, `ideas.md`, and the relevant observer reference notes reflect the landed direction.
- [ ] No regressions in `test/test_features.py`, `test/test_simulation_contracts.py`, `test/test_opencl_structs.py`, or `test/kernel_components/`.

## Key Refs

- `.design/reference/observer_readout_audit.md` — current readout contract plus the 2026-06 audit findings
- `.design/reference/observer_concept_audit.md` — deeper observer rationale and remaining design questions
- `.design/reference/compatibility_boundary_audit.md` — canonical vs compatibility surface policy after the audit
- `.design/reference/observer_solution_buffer_audit.md` — landed shared-history/helper contract
- `clode/observers/_definitions.py`, `clode/observers/types.py`, `clode/simulation/features.py` — Python observer surface and schema definitions
- `clode/kernels/observers/` — kernel-side family implementations
