# Single-Precision Numerics Note

Purpose: record the verified float32 failure modes, the current demonstration strategy, and the guardrails for future mitigation work.
Read when: touching compensated mean helpers, time-base updates, or public docs about single-precision accuracy.
Update when: the demonstration scripts change, helper adoption broadens, or the time representation changes.

## Verified points

- `runningMeanTime(...)` and the compensated integral path are algebraically related but not numerically equivalent in float32.
- The current public demonstration baseline is `examples/single_precision_accuracy.py`, which compares a naive incremental mean against compensated integral accumulation and compares direct time addition against Kahan-style and structured time updates.
- The public example stays NumPy-based on purpose so the formulas are easy to inspect and rerun; the actual OpenCL helper prototypes are pinned separately in `test/kernel_components/test_kernel_math.py`.
- The current stepper issue is representational, not stylistic. `ti += dt` and `ti = ti + dt` both perform the same float32 addition against the current absolute time.
- The live fixed-step path now reconstructs absolute time from `t0 + step * dt` with a step counter, while the adaptive steppers still derive the next absolute time from one float32 absolute-time value plus one float32 step size.
- Kahan-style or structured time updates improve long-horizon endpoint accuracy and reduce drift, but they do not manufacture representable intermediate times once `dt < ulp(t)`.
- Reconstructing time as `t = step * dt` or `t = t0 + step * dt` in float32 avoids repeated-add drift, but it is still subject to float32 spacing and, if the step counter is cast to float32, it loses unit-step exactness beyond `2^24`.
- The live fixed-step path now uses a 64-bit step counter, so loop-budget overflow is much less constraining than float32 time spacing for practical workloads.
- The new helper-layer prototypes split the time problem more cleanly: fixed-step kernels can use `t0 + step * dt` from a 64-bit step counter, while adaptive kernels need a two-float compensated time pair because there is no fixed-step counter shortcut.
- The compensated integral helper used by `basic` and `basicall` has a clear empirical justification from the public example cases; broader observer adoption should be justified the same way rather than assumed.
- Large-origin elapsed-time differences are a separate live failure mode for feature extraction: if the mean denominator is formed as `ti - t_start` in float32, the subtraction can be badly biased or collapse to zero even when the numerator is accumulated carefully.
- The current OpenCL path should be treated as standard round-to-nearest-even arithmetic unless a future experiment proves otherwise. Stochastic rounding may reduce bias in some low-precision accumulations, but it does not remove float32 representability limits and would complicate reproducibility.
- The current public examples now show three observer-safeguard cases that are actually live today: `min_amp` suppressing tiny oscillations, `threshold_2` Schmitt-trigger-style hysteresis suppressing ripple-driven chatter, and `dx_up_threshold` suppressing noisy shallow crossings in a slope-gated case.
- The current public example also shows that threshold-crossing timestamps have a clear alternative hierarchy on smooth coarse crossings: sampled endpoint times are crude, inverse-linear interpolation is already much better, and slope-aware Hermite interpolation can be dramatically more accurate when the crossing is monotone.
- A follow-up ripple-heavy crossing check shows why Hermite should stay prototype-only for now: when a coarse timestep contains multiple threshold crossings, the interpolation target is ambiguous and inverse-linear interpolation can be more robust than a higher-order single-crossing model.
- The current public example also includes a three-sample local-extremum prototype showing that quadratic-vertex fitting on the existing buffer geometry can materially outperform sample-argmax localization on coarse data.
- `clode/kernels/clODE_utilities.cl` now includes tested prototypes for Kahan-style time accumulation, fixed-step counter time reconstruction, inverse-linear threshold timestamps, slope-aware Hermite threshold timestamps, and bounded three-sample max/min helpers. The live fixed-step steppers now use the counter-reconstructed time path, `threshold_2` now uses inverse-linear timestamps for stored up/down threshold transitions, and `local_max` now uses the bounded three-sample extremum helpers, but broader helper adoption is still selective and incomplete.
- `min_imi` and `eps_dx` are still exposed on `ObserverParams`, but they are not yet wired into the current built-in observer kernels strongly enough to recommend as primary safeguards in public docs.

## Public-docs boundary

Public docs should be explicit about three things:

- where compensated mean accumulation demonstrably helps
- where feature-window origin matters because `ti - t_start` can fail independently of the compensated numerator
- where time compensation or structured time helps only with drift and endpoint accuracy
- where current observer timestamps are sample-based and what inverse-linear or slope-aware interpolation alternatives buy on smooth crossings
- where a stronger fix likely needs a richer time representation than one float32 absolute-time value
- public package docs should describe current supported behavior and compare it to alternative algorithms, not to the package's development history

Public docs should recommend starting autonomous feature windows near zero when that does not change the model semantics.

Avoid public wording that implies notation cleanup alone mitigates the time-base problem.

## Preferred implementation direction

- For future single-precision GPU paths, prefer a solver-owned dual-realtype compensated time representation by default rather than another single-float absolute-time variant.
- For fixed-step steppers, keep using a wide step counter plus `t0 + step * dt` reconstruction instead of repeated float32 absolute-time addition.
- Use the same richer time representation, or a solver-owned relative elapsed channel derived from it, for observer elapsed-time bookkeeping so time-weighted means do not depend on subtracting two large float32 absolute times.
- Treat broader compensated-mean rollout and richer time-base work as one coherent bookkeeping problem rather than as isolated helper swaps.
- Investigate selective adoption of the new shared three-sample local-extremum helper and threshold-timestamp interpolation utilities before changing each observer independently.
- None of the current helper candidates needs a general nonlinear system solver. The Hermite threshold prototype only uses a scalar Newton iteration on a cubic. Any future helper that depends on solving a coupled nonlinear system should be deferred with the later implicit-method work.

## Follow-on questions

- Which remaining observers actually benefit enough from compensated integral accumulation to justify the extra state fields?
- For non-autonomous problems, what is the lightest implementation that preserves the preferred direction above: dual-realtype solver time everywhere, a dual-realtype plus relative-elapsed split, or some narrower per-work-item `t0` plus relative-time model?
- Which of those options is worth carrying into the live kernels before any broader public API redesign?
- Beyond the current live 64-bit fixed-step counter path, is any more elaborate coarse/fine counter scheme actually worth the complexity?
- For threshold-based observers, is inverse-linear interpolation the right default first step, and where would a Hermite-style slope-aware helper remain robust enough to justify its extra complexity?
- Is there any practical OpenCL path to experiment with stochastic rounding while preserving clODE's reproducibility expectations, or should rounding-mode work stay out of scope?
- Which live kernels need more scrutiny for cancellation, threshold jitter near zero, or divisions by tiny ranges beyond the current mean and time demos?

## Other live hot spots

- Event detectors in `clode/kernels/observers/observer_local_maximum.clh`, `observer_threshold_2.clh`, `observer_neighborhood_1.clh`, and `observer_neighborhood_2.clh` rely on sign changes such as `dx > 0` then `dx < 0`; those comparisons are sensitive to slope noise near zero in float32.
- `observer_basic.clh` and `observer_basic_allVar.clh` still derive elapsed time as `ti - t_start`, so feature windows started at large absolute times remain a real precision hotspot even after numerator compensation.
- `observer_threshold_2.clh` now uses inverse-linear timestamps for stored up/down threshold transitions in the live kernel, using the later active boundary when both `x` and `dx` gates participate. Slope-aware Hermite interpolation should still stay prototype-only until its robustness on noisy crossings is better characterized.
- `observer_local_maximum.clh` now uses the shared bounded three-sample max/min helpers, but extrema on coarse or noisy traces remain a hotspot worth watching.
- Neighborhood observers normalize by trajectory ranges in `observer_neighborhood_1.clh` and `observer_neighborhood_2.clh`; tiny or poorly resolved ranges can amplify noise or push the computation toward division by a very small number.
- Interpolation helpers in `clode/kernels/clODE_utilities.cl` divide by `t1 - t0` and related expressions, so closely spaced or quantized time buffers are a natural place to watch for loss of significance.
- Any future variance, covariance, or higher-moment features should avoid naive subtract-large-sums formulas; the same catastrophic-cancellation pattern that motivated compensated means will matter there too.
