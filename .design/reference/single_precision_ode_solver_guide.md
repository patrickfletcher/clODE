# Single-Precision Floating-Point Arithmetic for Sensitive ODE Solvers

Purpose: float32 numerics guidance for sensitive solver and kernel work.
Read when: editing `clode/kernels/*`, solver numerics, reductions, tolerances, scaling, or mixed-precision behavior.
Update when: float32 guardrails, kernel numerics guidance, or recommended implementation patterns change.

## IEEE-754 `float32` essentials

- A `float32` has 1 sign bit, 8 exponent bits, and 23 explicit mantissa bits, or about 24 bits of effective precision with the hidden leading 1.
- Key quantities:

| Quantity | Approximate Value |
| --- | --- |
| Machine epsilon (`eps`) | `1.1920929e-7` |
| Decimal precision | ~7 significant digits |
| Largest finite value | `3.4028235e38` |
| Smallest normal value | `1.1754944e-38` |
| Smallest subnormal value | `1.4012985e-45` |

- `eps` is spacing near `1.0`; relative spacing scales with magnitude.
- Floating-point arithmetic is not associative:

```text
(a + b) + c != a + (b + c)
(a * b) * c != a * (b * c)
```

- It is often not distributive:

```text
a * (b + c) != a*b + a*c
```

- Do not assume algebraic rearrangements preserve numerics.

## Core failure modes

### 1. Catastrophic cancellation

- Subtracting nearly equal numbers destroys leading digits:

```text
x = 1.234567
y = 1.234566
x - y = 1e-6
```

- Important in finite differences, residual computations, Jacobian estimation, conservation-law corrections, and variance-like formulas.
- Mitigations: avoid near-equal subtraction; reformulate analytically; use compensated formulas; maintain numerically stable coordinates; prefer fused expressions over staged differences.
- Example:

```text
Bad:    sqrt(x+1) - 1
Better: x / (sqrt(x+1) + 1)
```

### 2. Accumulation error

- Repeated summation accumulates rounding error.
- Naive summation error grows roughly like `O(N * eps)`, or worse when magnitudes vary strongly.
- Critical in RHS accumulation, reduction kernels, norm computations, residuals, statistics, and long trajectories.
- Mitigations:

#### Kahan summation

- Track lost low-order bits:

```text
sum = 0
c = 0

for x in values:
        y = x - c
        t = sum + y
        c = (t - sum) - y
        sum = t
```

- Often worth the cost in sensitive kernels.

#### Pairwise reduction

- Prefer tree reductions over sequential accumulation; this is good for SIMD and GPU implementations.

#### Sum small terms first

- Ascending-magnitude summation reduces error.

### 3. Dynamic range problems

- Overflow: `exp(100)` may overflow in `float32`.
- Underflow: very small values collapse toward zero.
- Subnormals preserve gradual underflow but are slow on many architectures and are often flushed to zero (FTZ).
- Mitigations: nondimensionalize equations; rescale states and parameters; work in log space where appropriate; avoid unnecessary exponentials and powers; detect dangerous intermediate expressions.

### 4. Loss of significance in time integration

- ODE solvers repeatedly add small increments to larger states:

```text
y_new = y + dt * f
```

- If `|dt * f| << |y|`, the increment may partially or completely vanish.
- Major in long integrations and multiscale systems.
- Mitigations: scale variables to `O(1)`; use adaptive timesteps; avoid excessively small `dt`; monitor increment and stagnation ratios; consider mixed-precision accumulators.
- Useful diagnostic:

```text
if y_new == y:
        update vanished
```

### 5. Stiffness amplifies precision problems

- Stiff systems generate widely separated timescales, ill-conditioned Jacobians, cancellation in implicit solves, and large transient amplification.
- Single precision reduces stable operating margins.
- Consequences: false convergence, Newton stagnation, unstable linear solves, inaccurate eigenstructure.
- Mitigations: aggressive scaling; analytically derived Jacobians; preconditioning; condition monitoring; residual-based stopping criteria; avoid over-solving linear systems beyond float32 resolution.
- Do not set nonlinear tolerances below meaningful float32 accuracy.
- Typical lower practical bound: `~1e-6 to 1e-7`, depending on conditioning.

## ODE-specific numerical guidance

### Error tolerances

- Setting tolerances below float32 resolution is meaningless.
- Typical practical ranges:

| Quantity | Typical Range |
| --- | --- |
| Relative tolerance | `1e-4` to `1e-6` |
| Absolute tolerance | problem-scaled |

- Avoid `rtol = 1e-10` in pure float32.

### State scaling

- High-value intervention.
- Aim for state variables roughly near `O(1)`.
- Avoid mixed scales like `1e9 and 1e-12` inside the same state vector.
- Scaling improves conditioning, timestep selection, Jacobians, residuals, and nonlinear solves.

### Jacobians

- Finite-difference Jacobians are fragile in float32.
- Bad step sizes produce cancellation and pure noise derivatives.
- Prefer analytic Jacobians, automatic differentiation, or complex-step differentiation when applicable.
- Finite-difference perturbations must scale carefully.

### Linear algebra

- Ill-conditioned linear systems destroy float32 accuracy quickly.
- Watch for nearly singular matrices, badly scaled pivots, and subtractive elimination.
- Mitigations: scaling or preconditioning, pivoting, iterative refinement, mixed-precision solves.
- Condition numbers near `1/eps ≈ 1e7` already threaten meaningful accuracy.

## Mixed precision strategy

- Mixed precision is often the best compromise.
- Common pattern:

| Operation | Precision |
| --- | --- |
| state storage | float32 |
| vectorized RHS | float32 |
| reductions/norms | float64 |
| timestep controller | float64 |
| accumulators | float64 |
| linear solves | mixed |

- A few float64 accumulators can dramatically improve robustness.

## GPU / SIMD considerations

- Many accelerators fuse multiply-add (`FMA`), reorder operations, use non-deterministic reductions, and flush subnormals to zero.
- Results vary across compiler versions, optimization levels, and hardware targets.
- Do not assume bitwise reproducibility.
- If reproducibility matters, disable unsafe math optimizations, control reduction ordering, avoid atomics for sensitive reductions, and test across architectures.

## Compiler and optimization hazards

- Flags like `-ffast-math` may violate IEEE assumptions.
- Possible effects: reassociation, dropped NaN checks, reciprocal approximations, altered transcendental accuracy.
- Use these flags carefully.

## NaN and Inf hygiene

- Numerical failures often first appear as NaN, Inf, denormal explosions, or silent divergence.
- Add `isfinite(x)` checks in development builds.
- Monitor norms, timestep collapse, residual growth, and conservation-law drift.

## Practical design principles

1. Scale everything. Poor scaling is a dominant cause of float32 instability.
2. Avoid numerically hostile algebra. Equivalent mathematics does not imply equivalent numerics.
3. Use stable reductions. Never trust naive summation in sensitive kernels.
4. Do not over-resolve. Float32 cannot support arbitrarily tight tolerances.
5. Prefer robustness over nominal order. A theoretically higher-order method can lose to a lower-order but numerically stable implementation.
6. Monitor conditioning. Precision problems are often conditioning problems in disguise.

## High-value diagnostics

| Diagnostic | Purpose |
| --- | --- |
| `y_new == y` | vanished updates |
| residual stagnation | precision floor |
| timestep collapse | stiffness / instability |
| conservation drift | accumulated error |
| Jacobian condition estimate | solve reliability |
| NaN/Inf detection | hard failure |
| norm growth | instability |

## Common false assumptions

| False | Reality |
| --- | --- |
| "Float32 gives 7 correct digits everywhere." | Precision is relative and magnitude dependent. |
| "If the method is mathematically stable, the implementation is stable." | Floating-point effects can destabilize stable algorithms. |
| "Smaller timestep always improves accuracy." | Eventually roundoff dominates truncation error. There is often an optimal timestep. |
| "Two algebraically equivalent formulas behave the same numerically." | Expression structure strongly affects rounding behavior. |

## Recommended default practices

Default float32 practice:

- nondimensionalize equations
- scale states near `O(1)`
- use adaptive timesteps
- avoid finite-difference Jacobians when possible
- use compensated or pairwise summation
- use float64 for reductions and controllers
- monitor conditioning
- detect NaN and Inf aggressively
- test across hardware and compiler configurations
- validate against float64 reference solutions

## Important final perspective

- For many ODE systems, float32 failure is driven less by the nominal precision limit itself and more by poor scaling, ill-conditioning, unstable algebraic structure, reduction error, stiffness, and inappropriate tolerances.
- Well-designed float32 solvers can be remarkably robust.
- Poorly scaled or numerically careless implementations can fail catastrophically even in float64.
