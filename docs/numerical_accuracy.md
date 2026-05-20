# Numerical Accuracy in Single Precision

clODE runs efficiently in single precision, but float32 arithmetic has several predictable limits that matter for long integrations and online feature extraction:

- time-weighted mean updates can stop responding to late small contributions
- elapsed-time differences can lose significance when both endpoints are large
- repeated time updates can drift or stall once $dt$ is small relative to the spacing between adjacent representable times
- event detectors and local-extremum estimates can overreact to chatter or coarse sampling unless the observer configuration matches the signal

The reproducible script `examples/single_precision_accuracy.py` exercises the same float32 update patterns used by the current kernels and observers.

## The Core Float32 Facts

IEEE-754 float32 has 24 bits of significand precision, so the unit roundoff is

$$
u = 2^{-24} \approx 5.96 \times 10^{-8}.
$$

Near a value with magnitude about $2^e$, the gap between adjacent float32 numbers is approximately

$$
\operatorname{ulp}(x) \approx 2^{e-23}.
$$

That matters in two places.

First, a recurrence such as

$$
m_n = m_{n-1} + (x_n - m_{n-1}) \frac{\Delta t_n}{T_n}
$$

eventually applies updates that are smaller than the local float32 spacing of $m_n$.

Second, direct time updates such as

$$
t_{n+1} = \operatorname{fl}(t_n + \Delta t)
$$

become unreliable once $\Delta t$ is small relative to $\operatorname{ulp}(t_n)$. The notation does not matter here: `ti += dt` and `ti = ti + dt` perform the same floating-point addition and rounding.

## Compensated Mean Accumulation

The `basic` and `basicall` observers accumulate

$$
I_n = \sum_{k=1}^{n} \Delta t_k x_k
$$

with Neumaier compensation and compute the mean once at the end as $I_n / T_n$. A useful comparison baseline is the naive incremental mean recurrence `runningMeanTime(...)`.

The example script compares those two formulas on two synthetic float32 workloads:

- a half-window step from $1.0$ to $1.0001$
- a long smooth oscillation, $1 + 10^{-4}\sin(t)$

With the default script parameters:

- the half-window step has exact mean $1.00005$
- `runningMeanTime` stays at exactly `1.0`, losing the late offset entirely
- the compensated integral path returns about `1.00005007`, keeping the error near $7 \times 10^{-8}$
- on the smooth oscillation, `runningMeanTime` drifts by about $7.1 \times 10^{-6}$ while the compensated integral path stays near $7.5 \times 10^{-8}$

This is the clearest current empirical justification for the compensated helper already used in `basic` and `basicall`.

## Feature Windows and Elapsed-Time Origin

Compensating the numerator is only part of the mean problem. The denominator still needs a trustworthy elapsed time.

The example script adds a large-origin stress test for a compensated mean where the numerator is accumulated in relative time, but the denominator is formed either as:

- `ti - t_start` in float32
- a separately accumulated relative elapsed time

With `dt = 0.01` over a 100-unit feature window:

- at `t0 = 0`, both paths return about `1.0000205`
- at `t0 = 1e4`, the `ti - t_start` path already biases the mean to about `1.02405` because the elapsed time snaps to `97.65625` instead of about `100`
- at `t0 = 1e6`, the `ti - t_start` path collapses completely because the float32 absolute time no longer advances relative to `t_start`

This is why feature extraction has an additional user-side rule in single precision: when the problem is autonomous and absolute time is only a label, prefer feature windows that start near `t = 0`. If the model is non-autonomous and absolute time matters semantically, treat the issue as a time-representation limitation rather than expecting numerator compensation alone to repair it.

## Time Accumulation: What Helps, What Does Not

The same example script compares three float32 time-update strategies:

- direct repeated addition
- Kahan-style compensated addition
- structured reconstruction, `t0 + step * dt`

### Zero-origin drift

Starting from `t0 = 0`, `dt = 0.01`, and 200,000 steps:

- direct repeated addition ends about `1.64` time units away from the float64 reference
- Kahan-style accumulation reaches the correct endpoint in this test
- `t0 + step * dt` also reaches the correct endpoint in this test

This shows that even when $dt$ is representable, naive float32 accumulation can drift badly over long runs.

### Large-origin stall

Starting from `t0 = 1e6`, `dt = 0.01`, and 10,000 steps:

- `ulp(t0)` is `0.0625`, which is already larger than `dt`
- direct repeated addition stalls on the first step and never advances
- Kahan-style and structured updates recover the final endpoint, but only by moving in `0.0625` jumps

Compensation and structured reconstruction improve endpoint accuracy, but they do not create representable times between adjacent float32 numbers. Once $dt < \operatorname{ulp}(t)$, exact per-step absolute-time resolution is unavailable in a single float32 time variable.

### Long zero-origin horizon

The same script also pushes the zero-origin case to 30,000,000 steps.

- direct repeated addition first produces a zero increment at about step `22918664` and finishes at `262144.0` instead of `300000.0`
- Kahan-style and structured updates both first show zero increments at about step `13107201`, but they still reach the correct endpoint in this test
- this means Kahan-style and structured time control drift much better than naive accumulation, but they still quantize into larger jumps once $dt$ falls below the local float32 spacing
- if the step counter itself is cast to float32, unit-step exactness is also lost beyond $2^{24} \approx 1.68 \times 10^7$

So `t = t0 + step * dt` is a bounded-error float32 strategy, not a complete fix for absolute-time representability.

For fixed-step steppers, this is now the live current strategy: keep a 64-bit integer step counter and reconstruct absolute time from `t0 + step * dt` instead of from repeated float32 addition. Adaptive steppers do not have that shortcut because `dt` changes from step to step, so their stronger prototype direction is still a dual-float compensated time pair rather than another single-float absolute-time variant.

## How To Run The Demonstration

Run the example from a source checkout:

```python
python examples/single_precision_accuracy.py
```

If `matplotlib` is available, the script also produces a summary figure.

## Observer Safeguards for Oscillatory Features

Some observer modes are intentionally tuned for oscillatory behavior. In those workflows, the right response to numerical fragility is often not only a different arithmetic scheme, but also a user-configurable guard that says which events are meaningful enough to count.

### Amplitude floor with `min_amp`

On a tiny sinusoid with total range about $10^{-3}$:

- a `threshold_2`-style event count with `min_amp = 0` reports 20 up-events across the sampled window
- the same trace with `min_amp = 2 \times 10^{-3}` reports 0 events

This encodes the smallest oscillation that is physically or scientifically meaningful for the workflow.

### Schmitt-trigger-style hysteresis with separate up/down thresholds

On a coarse sine wave with a high-frequency ripple superimposed:

- a single threshold (`x_up = x_down = 0.5` in normalized coordinates) reports 100 up-events
- a clear hysteresis gap (`x_up = 0.65`, `x_down = 0.35`) returns the expected coarse-cycle count of 20 up-events

Using `x_up > x_down` is intentional here. It is the same basic idea as a Schmitt trigger: do not treat every ripple near the boundary as a new event.

### Derivative thresholds as a slope gate

On a composite waveform where value thresholds alone still overcount crossings:

- a `threshold_2`-style count with `dx_up = 0` reports 100 up-events
- adding `dx_up = 0.9` reduces that to 20 up-events

Derivative thresholds are therefore a real current safeguard for some noisy traces, especially when value hysteresis alone is not selective enough. They are still signal-dependent, so validate them on representative data before launching a large sweep.

## Threshold-Crossing Timestamp Alternatives

The current `threshold_2` kernel records transition times at the sampled endpoint of the crossing step. The example script compares that current sample-time choice against two interpolation alternatives on a coarse smooth crossing:

- sampled crossing time
- inverse-linear interpolation between the two bracketing samples
- cubic Hermite interpolation using the same endpoint values and endpoint slopes

On an upward threshold crossing of a coarsely sampled sine wave:

- the sampled timestamp has mean error about `9.61e-2`
- inverse-linear interpolation reduces that to about `1.97e-3`
- slope-aware Hermite interpolation reduces it further to about `1.33e-6`

This is enough evidence to justify an explicit threshold-crossing helper in the shared kernel math layer. It also shows a useful distinction: inverse-linear interpolation is the low-risk option when you only want a better timestamp, while slope-aware Hermite interpolation is an attractive prototype when the crossing is smooth and monotone but still needs stronger robustness checks before broad use in noisy event detectors.

That caveat is real rather than theoretical. On ripple-heavy traces, one coarse timestep can contain multiple threshold crossings, which makes any single-crossing interpolation model ambiguous. That is why inverse-linear interpolation is the conservative first live candidate and slope-aware Hermite interpolation remains prototype-only for now.

## Local-Extremum Localization from a Three-Sample Buffer

The current `local_max` observer localizes extrema by picking the sampled maximum from its three-sample buffer. The example script also includes a prototype based on that same geometry.

On a coarsely sampled sine wave:

- taking the buffered sample maximum directly gives a mean peak-time error of about `5.15e-2` and a mean peak-value error of about `1.78e-3`
- fitting a quadratic vertex on the same three samples reduces the mean peak-time error to about `3.57e-4` and the mean peak-value error to about `2.66e-5`

That does not change current package behavior by itself, but it does show why local-extremum detection remains a live accuracy hotspot. If peak timing matters today, the safe user-side choices are a smaller `dt` or double precision.

## User Recommendations

- The compensated integral helper used by `basic` and `basicall` has a clear empirical justification.
- For autonomous feature extraction, prefer windows that start near `t = 0` when possible; large absolute times can break float32 elapsed-time differences even when the numerator is accumulated carefully.
- For non-autonomous systems where absolute time matters, prefer double precision today or shorter windows that keep `t` near the scale you need.
- The stepper time-base issue is a representability problem, not a notation problem.
- For fixed-step steppers, the live step-counter time base is the cleanest single-precision alternative to repeated addition, but it is still a float32 strategy and therefore still quantized by float32 spacing.
- For adaptive steppers, a stronger fix needs a dual-float compensated time representation or double precision because there is no fixed-step counter shortcut.
- For `threshold_2`, use a real hysteresis gap when crossings chatter, and consider derivative thresholds when noisy shallow crossings still slip through.
- Treat current threshold timestamps as sampled times; when coarse timestamp accuracy matters, inverse-linear interpolation is the first alternative worth evaluating.
- Set `min_amp` above the numerical or measurement floor you want to ignore and below the smallest oscillation you still care about.
- For `local_max`, use a smaller `dt` or double precision when accurate peak timing matters.
- Validate event-sensitive settings on one representative trajectory before launching a large ensemble sweep.

For the runnable script inventory, see [examples.md](examples.md). For throughput-oriented timing notes, see [performance_notes.md](performance_notes.md).
