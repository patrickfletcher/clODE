# Fixed-Step Endpoint Bug Audit

## Scope

This note re-audits the fixed-step endpoint issue after changing the main time-loop guards from `<=` to `<`.

The question here is narrower than in the first draft: with fixed-step methods retaining constant step size, why do `tf` and `xf` still sometimes end up beyond `tspan[1]`?

Files inspected for this revision:

- `clode/cpp/transient.cl`
- `clode/cpp/trajectory.cl`
- `clode/cpp/features.cl`
- `clode/cpp/initializeObserver.cl`
- `clode/cpp/steppers/fixed_explicit_step.clh`
- `clode/cpp/steppers/fixed_explicit_Euler.clh`
- `clode/cpp/steppers/fixed_explicit_RK4.clh`
- `clode/cpp/steppers/fixed_explicit_Trapezoidal.clh`
- `clode/cpp/steppers/adaptive_explicit_step.clh`
- `clode/cpp/observers/observer_basic_allVar.clh`
- `clode/cpp/CLODE.cpp`

## Executive Summary

The current fixed-step endpoint issue is best described as a combination of:

1. repeated floating-point accumulation in `ti = ti + dt`
2. a pre-step loop guard of the form `while (ti < tspan[1])`
3. reporting the post-step state and time after the loop exits

After the guard change from `<=` to `<`, the bug is no longer explained by an unconditional extra step at an exact endpoint. For exact binary cases such as `dt = 0.25` on `(0, 1)`, the run now stops cleanly at `tf = 1.0`.

The remaining overshoot happens when accumulated rounding leaves `ti` just below `tspan[1]` at the start of an iteration. The loop then correctly takes one more fixed step, and the returned `xf` and `tf` correspond to the state after that step. That is why the issue persists with `<`.

This is not a reason to shrink the last step. For a fixed-step solver, the step size should remain fixed. Any future fix should preserve constant-step semantics.

## Current Mechanism

### 1. The fixed-step loop checks the time before the step

The relevant kernels still follow the same pattern:

- test `ti < tspan[1]`
- take one full fixed step
- store the updated `ti` and `xi`

That means the condition only governs whether another step starts. It does not guarantee that the post-step time remains inside the requested interval.

### 2. Fixed-step methods always advance by the full `dt`

This is the correct design principle for fixed-step methods.

The wrapper in `fixed_explicit_step.clh` calls the stepper with `*dt`, and the fixed-step methods update time by that same amount:

- Euler: `*ti += dt`
- RK4: `*ti = *ti + dt`
- Trapezoidal/Heun: `*ti = *ti + dt`

There should be no adaptive-style end-step shrink in this path.

### 3. Floating-point accumulation decides whether one more step starts

When `dt` is not exactly representable in the chosen `realtype`, repeated additions can leave the running time slightly below the requested endpoint even when the mathematically ideal time would be exactly on it.

If the loop sees:

```text
ti = 0.9999994
tspan[1] = 1.0
```

then `ti < tspan[1]` is still true, so one more full step is taken and the reported final time becomes approximately:

```text
tf = 1.0099994
```

That is the behavior currently observed.

## Device Measurements

Using platform 1 device 0 with the stable linear test problem, the current fixed-step RK4 endpoint behavior is:

| `dt` | reported `tf` | endpoint error `tf - 1.0` |
| --- | --- | --- |
| `0.25` | `1.0` | `0.0` |
| `0.1` | `1.0000001192092896` | `1.19e-7` |
| `0.025` | `1.0249996185302734` | `2.50e-2` |
| `0.01` | `1.009999394416809` | `9.999e-3` |

These measurements support three conclusions:

- the old claim that `<` is insufficient because the loop always steps once past the endpoint in exact arithmetic is no longer correct
- endpoint overshoot now depends on floating-point accumulation and representability of `dt`
- the endpoint error generally becomes small when `dt` is sufficiently small for the interval and the accumulated time remains close to the target

The `0.025` case shows that smaller `dt` does not guarantee a smaller endpoint error for every decimal value. What matters is the interaction between `dt`, the interval length, and the floating-point representation used in accumulation.

## What The Loop Is Actually Returning

This is the key point behind the persistent bug.

The loop condition is checked using the current time before a step starts. After the step, the solver keeps the new `ti` and `xi`. When the next loop check fails, those stored values are already the post-step endpoint.

So `xf` and `tf` are not "the last state inside the interval". They are "the state and time after the last accepted fixed step".

That explains why changing `<=` to `<` reduces the issue but does not eliminate it.

## Adaptive Contrast

The adaptive path is different for a good reason.

Adaptive methods already vary step size by design, so shrinking the final accepted step to respect the endpoint is normal there. The adaptive wrapper in `adaptive_explicit_step.clh` clamps the proposed step against the remaining interval.

That logic should remain confined to adaptive steppers. It is not an appropriate fix for fixed-step methods.

## Downstream Effects

### 1. Final-time reporting

`get_tf()` can report a value slightly beyond the requested endpoint. Tests that compare against exact solutions should therefore use the returned time when validating the returned state.

### 2. Trajectory sampling

Trajectory outputs store the time and state after each accepted step. If the final accepted step crosses the endpoint, the last stored sample can land beyond the requested window.

### 3. Feature summaries

The `basicall` observer updates from the actual accepted step sequence. If one more step is accepted because accumulated time is still below the endpoint, the observer statistics reflect that extra sample.

### 4. Continuation semantics

Continuation remains sensitive because `xf` and `tf` correspond to the post-step state/time, while utilities such as `shiftTspan()` are defined in terms of the requested span rather than the reported final time.

## Separate Observation: `basicall` Continuation Still Looks Distinct

The feature-continuation mismatch still appears to involve more than endpoint handling alone.

`observer_basic_allVar.clh` contains:

```c
static inline void finalizeObserverData(...)
{
    realtype T = *ti - tspan[0];
    od->t_start -= T;
}
```

That update remains suspicious for split-window continuation because it mutates the observer time base across runs. This looks like a separate issue from the fixed-step endpoint behavior.

## Practical Recommendations

### 1. Do not add last-step shrinking to fixed-step solvers

That would blur the distinction between fixed-step and adaptive methods.

### 2. If this is fixed later, preserve constant-step semantics explicitly

Two plausible directions are:

1. keep an integer step counter and reconstruct time from `t0 + step * dt`
2. keep accumulated time but use compensated summation such as Kahan summation

Either of those is more consistent with fixed-step semantics than adaptive-style end-step clipping.

### 3. Treat endpoint drift as a documented backend quirk for now

Since further solver redesign is likely anyway when adding fixed-step implicit and adaptive implicit methods, the most practical near-term move is to document the current behavior and defer a structural fix.

### 4. Write deterministic tests with sufficiently small fixed steps

When validating exact deterministic solutions, use a `dt` small enough that the ODE discretization error is clearly below the desired assertion tolerance. That keeps the tests focused on solver correctness rather than on intentionally crude coarse-step behavior.

## Bottom Line

With the current `<` guard, the fixed-step endpoint issue is no longer primarily a loop-boundary bug by itself. The remaining problem is that floating-point accumulation can leave `ti` slightly below `tspan[1]`, causing one more full fixed step to be taken, and the solver then reports the post-step `xf` and `tf`.

That is a real quirk of the current backend, but the correct future fix should retain true fixed-step semantics rather than borrowing adaptive end-step behavior.