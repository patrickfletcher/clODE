# Single-Precision Numerics Note

Purpose: record the verified float32 failure modes, the current demonstration strategy, and the guardrails for future mitigation work.
Read when: touching compensated mean helpers, time-base updates, or public docs about single-precision accuracy.
Update when: the demonstration scripts change, helper adoption broadens, or the time representation changes.

## Verified points

- `runningMeanTime(...)` and the compensated integral path are algebraically related but not numerically equivalent in float32.
- The current public demonstration baseline is `examples/single_precision_accuracy.py`, which compares a naive incremental mean against compensated integral accumulation and compares direct time addition against Kahan-style and structured time updates.
- The public example stays NumPy-based on purpose so the formulas are easy to inspect and rerun; the actual OpenCL helper prototypes are pinned separately in `test/kernel_components/test_kernel_math.py`.
- The current stepper issue is representational, not stylistic. `ti += dt` and `ti = ti + dt` both perform the same float32 addition against the current absolute time.
- The live fixed-step and adaptive paths both accumulate solve-relative time with a compensated two-realtype elapsed pair and reconstruct stage and endpoint times from `t_origin + elapsed`.
- The Kahan-style compensated time helpers now use internally consistent readback semantics: the stored correction term is subtracted back out when reconstructing elapsed time or absolute time.
- The live stepper wrappers now return the accepted step width directly to observer paths instead of asking `features.cl` to recover it as `elapsed_new - elapsed_old`; differencing large float32 elapsed totals can quantize small accepted steps away even when the compensated pair itself is still tracking the solve correctly.
- Fixed Heun and RK4 now build internal stage times from `t_origin + elapsed + fractional_dt` through the same compensated helper family used by adaptive steppers instead of adding fractions of `dt` to an already rounded absolute `ti`.
- Fixed-step `do_step` implementations now own compensated endpoint advancement and endpoint slope refresh, so precomputed `k1` means the current accepted-state slope on entry and the accepted-endpoint slope on return across the live stepper paths.
- Kahan-style or structured time updates improve long-horizon endpoint accuracy and reduce drift, but they do not manufacture representable intermediate times once `dt < ulp(t)`.
- Reconstructing time as `t = step * dt` or `t = t0 + step * dt` in float32 avoids repeated-add drift, but it is still subject to float32 spacing and, if the step counter is cast to float32, it loses unit-step exactness beyond `2^24`; that remains useful as a comparison baseline even though it is no longer the live kernel path.
- The current public comparison baseline keeps the step counter wide enough that loop-budget overflow is much less constraining than float32 time spacing for practical workloads.
- The helper-layer prototypes still support both comparison baselines cleanly, but the live kernels now use the same two-realtype compensated elapsed pair across fixed-step and adaptive steppers because the current direct evidence does not show a clear fixed-step advantage for keeping a separate structured path.
- The compensated integral helper used by `basic` and `basicall` has a clear empirical justification from the public example cases; broader observer adoption should be justified the same way rather than assumed.
- Large-origin elapsed-time differences are a separate live failure mode for feature extraction: if a mean denominator or period or duration is formed from subtracting large float32 absolutes, the result can be badly biased or collapse to zero even when the numerator is accumulated carefully. The live `basic`, `basicall`, `threshold_2`, `local_max`, `neighborhood_1`, and `neighborhood_2` paths now avoid that failure by tracking relative elapsed time directly.
- The current OpenCL path should be treated as standard round-to-nearest-even arithmetic unless a future experiment proves otherwise. Stochastic rounding may reduce bias in some low-precision accumulations, but it does not remove float32 representability limits and would complicate reproducibility.
- The current public examples now show three observer-safeguard cases that are actually live today: `min_amp` suppressing tiny oscillations, `threshold_2` Schmitt-trigger-style hysteresis suppressing ripple-driven chatter, and `dx_up_threshold` suppressing noisy shallow crossings in a slope-gated case.
- The current public example also shows that threshold-crossing timestamps have a clear alternative hierarchy on smooth coarse crossings: sampled endpoint times are crude, inverse-linear interpolation is already much better, and slope-aware Hermite interpolation can be dramatically more accurate when the crossing is monotone.
- A follow-up ripple-heavy crossing check shows why Hermite should stay prototype-only for now: when a coarse timestep contains multiple threshold crossings, the interpolation target is ambiguous and inverse-linear interpolation can be more robust than a higher-order single-crossing model.
- The current public example also includes a three-sample local-extremum prototype showing that quadratic-vertex fitting on the existing buffer geometry can materially outperform sample-argmax localization on coarse data.
- The live `neighborhood_return` and `nhood2` observers now keep their sampled anchors `x0` but refine stored exit times by linearly interpolating the full normalized state between the last inside sample and the first outside sample of the neighborhood ball.
- `clode/kernels/clODE_utilities.cl` now includes tested compensated-time helpers plus inverse-linear threshold timestamps, slope-aware Hermite threshold timestamps, and bounded three-sample max/min helpers. The public example still includes structured `t0 + step * dt` as a comparison baseline, `threshold_2` now uses inverse-linear timestamps for stored up/down threshold transitions, and `local_max` now uses the bounded three-sample extremum helpers, but broader helper adoption is still selective and incomplete.
- `min_imi` and `eps_dx` are still exposed on `ObserverParams`, but they are not yet wired into the current built-in observer kernels strongly enough to recommend as primary safeguards in public docs.

## Public-docs boundary

Public docs should be explicit about three things:

- where compensated mean accumulation demonstrably helps
- where feature-window origin still matters because float32 absolute timestamps remain quantized even after elapsed-time bookkeeping is separated from them
- where time compensation or structured time helps only with drift and endpoint accuracy
- where current observer timestamps are sample-based, where neighborhood-based timestamps now use sampled anchors plus full-state exit interpolation, and what inverse-linear or slope-aware interpolation alternatives buy on smooth crossings
- where a stronger fix likely needs a richer time representation than one float32 absolute-time value
- public package docs should describe current supported behavior and compare it to alternative algorithms, not to the package's development history

Public docs should recommend starting autonomous feature windows near zero when that does not change the model semantics.

Avoid public wording that implies notation cleanup alone mitigates the time-base problem.

## Preferred implementation direction

- Keep the solver-owned dual-realtype compensated time representation that is now live across fixed-step and adaptive steppers rather than reintroducing another single-float or structured fixed-step variant; the current direct demos show no clear fixed-step advantage that justifies two live time-base stories.
- Keep using solver-owned relative elapsed bookkeeping derived from the richer time base for observer elapsed-time paths so time-weighted means, periods, and durations do not depend on subtracting two large float32 absolute times.
- Keep fixed and adaptive stepper ownership aligned: wrappers should stay minimal, while each concrete `do_step` owns the stage geometry, compensated endpoint advance, and endpoint slope handoff that define the method.
- Keep treating accepted step width and next-step proposal as distinct quantities in adaptive wrappers. The former is observer input; the latter is controller state.
- Keep preserving stepper failure status at the solver boundary even before the fuller solver-state model lands; observers should ingest accepted step widths and timestamps, not own failure policy.
- Treat broader compensated-mean rollout and richer time-base work as one coherent bookkeeping problem rather than as isolated helper swaps.
- Investigate selective adoption of the new shared three-sample local-extremum helper and threshold-timestamp interpolation utilities before changing each observer independently.
- None of the current helper candidates needs a general nonlinear system solver. The Hermite threshold prototype only uses a scalar Newton iteration on a cubic. Any future helper that depends on solving a coupled nonlinear system should be deferred with the later implicit-method work.

## Follow-on questions

- Which remaining observers or outputs benefit enough from richer relative-time bookkeeping to justify additional state fields?
- How much of the heavy-observer cost is true semantic state, such as retained event timestamps and per-variable buffers, versus avoidable private scratch or duplicated sample-time geometry?
- For non-autonomous problems, when is the current solve-local compensated pair enough, and when would a persisted per-work-item current-time or `t0` model pay for itself?
- Which of those options is worth carrying into the live kernels before any broader public API redesign?
- Beyond the current live 64-bit fixed-step counter path, is any more elaborate coarse/fine counter scheme actually worth the complexity?
- For threshold-based observers, is inverse-linear interpolation the right default first step, and where would a Hermite-style slope-aware helper remain robust enough to justify its extra complexity?
- Is there any practical OpenCL path to experiment with stochastic rounding while preserving clODE's reproducibility expectations, or should rounding-mode work stay out of scope?
- Which live kernels need more scrutiny for cancellation, threshold jitter near zero, or divisions by tiny ranges beyond the current mean and time demos?

## Other live hot spots

- Event detectors in `clode/kernels/observers/observer_local_maximum.clh`, `observer_threshold_2.clh`, `observer_neighborhood_1.clh`, and `observer_neighborhood_2.clh` rely on sign changes such as `dx > 0` then `dx < 0`; those comparisons are sensitive to slope noise near zero in float32.
- In the heavier event detectors, retained event arrays and per-variable buffers are usually a much larger per-work-item cost than the added `elapsedbuffer[3]` fields themselves, so `max_event_timestamps` and observer choice matter more than a single extra three-sample time buffer.
- Stored absolute timestamps in `observer_threshold_2.clh`, `observer_local_maximum.clh`, `observer_neighborhood_1.clh`, and `observer_neighborhood_2.clh` remain limited by float32 spacing at large absolute times even though the corresponding periods or durations now use relative elapsed bookkeeping.
- `observer_threshold_2.clh` now uses inverse-linear timestamps for stored up/down threshold transitions in the live kernel, using the later active boundary when both `x` and `dx` gates participate. Slope-aware Hermite interpolation should still stay prototype-only until its robustness on noisy crossings is better characterized.
- `observer_local_maximum.clh` now uses the shared bounded three-sample max/min helpers, but extrema on coarse or noisy traces remain a hotspot worth watching.
- Neighborhood observers normalize by trajectory ranges in `observer_neighborhood_1.clh` and `observer_neighborhood_2.clh`; tiny or poorly resolved ranges can amplify noise or push the computation toward division by a very small number.
- `observer_neighborhood_return.clh` and `observer_neighborhood_2.clh` now refine normalized-ball exit times from the last inside and first outside samples, but their anchors `x0` still come from the first sampled point below the anchor threshold, so coarse anchor steps remain a live accuracy limit.
- Interpolation helpers in `clode/kernels/clODE_utilities.cl` divide by `t1 - t0` and related expressions, so closely spaced or quantized time buffers are a natural place to watch for loss of significance.
- A leaner observer-time architecture may still be possible for `threshold_2` and `local_max`, for example by carrying elapsed-only or local-`dt` sample geometry and reconstructing absolute timestamps only when needed, but that should stay a measured design exercise because the current dual absolute/elapsed bookkeeping is what prevents large-origin subtraction from corrupting periods and durations.
- Any future variance, covariance, or higher-moment features should avoid naive subtract-large-sums formulas; the same catastrophic-cancellation pattern that motivated compensated means will matter there too.
