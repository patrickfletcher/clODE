# ODE Event Interpolation and Refinement

Purpose: capture practical interpolation and event-refinement guidance for observer and trajectory workflows, with current emphasis on linear and three-point quadratic methods.
Read when: choosing event-time or event-state refinement methods, deciding minimum sample buffers for observer families, or updating shared interpolation helpers.
Update when: clODE adopts new interpolation helpers or event-detection semantics change materially.

## Fast path

- For current observer implementation scope, this note is guidance input to `.design/reference/observer_solution_buffer_audit.md`.
- Current practical stance:
  - inverse-linear interpolation remains the default threshold-crossing timestamp refinement baseline.
  - bounded three-point quadratic interpolation is the preferred higher-order follow-through for threshold crossing and extremum refinement on accepted-step buffers.
- Cubic Hermite interpolation should not be treated as current follow-on direction for observer event refinement.

## Bottom line

- Threshold crossings and extrema should be treated differently:
  - threshold crossings: start with inverse-linear interpolation, then evaluate three-point quadratic refinement for better error order where buffer geometry supports it.
  - extrema: continue using bounded three-point quadratic refinement on accepted-step buffers.
- The currently preferred refinement progression is:
  1. sampled event detection semantics
  2. inverse-linear threshold timestamp refinement (already deployed in current kernels where applicable)
  3. three-point quadratic threshold time/state refinement using accepted-step buffers
- Any higher-complexity interpolants should require strong empirical evidence beyond these two stages.

## Buffer geometry guidance

- `K=2` accepted-step buffers are sufficient for sampled detection and inverse-linear threshold time refinement.
- `K=3` accepted-step buffers are required for three-point quadratic refinement.
- Where a family already carries `K=3` state/slope history, quadratic refinement can often be adopted without introducing new rolling-buffer shape.

## Design constraints

- Keep detection semantics unchanged when adding refinement.
- Keep refinement local to event output computation; do not silently shift trigger semantics.
- Prefer bounded interpolation helpers with explicit edge-case handling.
- Require focused tests for coarse-step behavior, near-flat slopes, and multi-crossing intervals.

## Historical note

Earlier Hermite-forward exploration is archived for reference only:

- `.design/archived/interpolation_policy_refresh_2026_05_31/ode_event_interpolation_note.md`

## Evidence anchors

- `clode/kernels/clODE_utilities.cl`
- `clode/kernels/observers/observer_threshold_crossing.clh`
- `clode/kernels/observers/observer_normalized_schmitt_trigger.clh`
- `clode/kernels/observers/observer_local_maximum.clh`
- `clode/kernels/observers/observer_normalized_neighborhood_return.clh`
- `test/kernel_components/test_kernel_math.py`
- `test/test_features.py`
