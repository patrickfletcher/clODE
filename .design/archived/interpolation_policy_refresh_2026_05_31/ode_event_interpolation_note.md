# ODE Event Interpolation and Refinement

Purpose: capture practical interpolation and event-refinement guidance for observer and trajectory workflows, including tradeoffs among linear, quadratic, and Hermite approaches.
Read when: choosing event-time or event-state refinement methods, deciding minimum sample buffers for observer families, or evaluating whether higher-order interpolation is worth the complexity.
Update when: clODE adopts new interpolation helpers, event-detection semantics change materially, or dense-output strategy changes.

## Fast path

Stop after this section unless you are evaluating interpolation formulas in detail.

- For current observer implementation scope, this note is guidance input to `.design/reference/observer_solution_buffer_audit.md` rather than a standalone implementation plan.
- Default practical stance: linear/inverse-linear threshold timing and bounded three-sample extrema stay baseline; cubic Hermite is a targeted upgrade path where endpoint slopes are reliable.
- This note does not by itself require broad kernel refactors or immediate higher-order dense-output rollout.

## Bottom line

- For smooth ODE solutions, **cubic Hermite interpolation** is usually the best practical choice.
- If you store:
  - time `t`
  - state `x`
  - derivative `dx/dt`

  at the **two most recent points**, you already have enough information for accurate interpolation and event refinement.
- Threshold crossings are generally easier and more accurate to localize than extrema.
- Extrema are intrinsically ill-conditioned when the derivative becomes very small over a broad interval.
- Linear interpolation is robust but low accuracy.
- Quadratic interpolation using 3 state samples is simple and good for extrema detection.
- Higher-order interpolation is usually unnecessary unless the ODE solver itself is high-order and already provides dense output coefficients.
- Many production ODE solvers internally use dense-output interpolants closely related to cubic Hermite interpolation.

---

# Problem setup

Suppose the ODE solver produces:

\[
\dot{x} = f(x,t)
\]

and stores a rolling buffer containing recent values of:

- \(t_i\)
- \(x_i\)
- \(\dot{x}_i\)

The goal is to:

1. Detect events:
   - threshold crossing
   - local extremum

2. Refine the event:
   - event time
   - event state

---

# 1. Linear interpolation

## Buffer requirements

### Threshold crossing
- Detection:
  - 2 state samples
  - \((t_0,x_0)\), \((t_1,x_1)\)
- Refinement:
  - same 2 samples

### Extremum
- Detection:
  - 3 state samples typically needed
- Refinement:
  - poor quality

---

## Method

Assume:

\[
x(t) \approx x_0 + \frac{x_1-x_0}{h}(t-t_0)
\]

where:

\[
h=t_1-t_0
\]

Threshold crossing estimate:

\[
t^* = t_0 + h\frac{x^*-x_0}{x_1-x_0}
\]

---

## Pros

- Extremely simple
- Cheap
- Robust
- Preserves monotonicity

---

## Cons

- Only first-order accurate
- Poor for extrema
- Ignores derivative information
- Often wastes the accuracy of higher-order ODE solvers

---

# 2. Quadratic interpolation (3 points)

## Buffer requirements

### Threshold crossing
- Detection:
  - 2 points sufficient
- Refinement:
  - 3 state samples improves accuracy

### Extremum
- Detection:
  - 3 state samples required
- Refinement:
  - same 3 samples

---

## Method

Fit a parabola through:

\[
(t_{-1},x_{-1}),\ (t_0,x_0),\ (t_1,x_1)
\]

For equally spaced samples, extremum time estimate:

\[
t^* \approx t_0
+
\frac{h}{2}
\frac{x_{-1}-x_1}
{x_{-1}-2x_0+x_1}
\]

---

## Pros

- Simple
- Better extremum localization
- Naturally matches local curvature near extrema
- No derivative storage needed

---

## Cons

- No derivative continuity
- Less accurate than Hermite interpolation
- More sensitive to uneven timestep spacing

---

# 3. Cubic Hermite interpolation (recommended)

## Buffer requirements

### Threshold crossing
- Detection:
  - 2 points
  - sign change in \(x-x^*\)
- Refinement:
  - 2 points with derivatives:
    - \((t_0,x_0,\dot{x}_0)\)
    - \((t_1,x_1,\dot{x}_1)\)

### Extremum
- Detection:
  - 2 derivative samples sufficient
  - sign change in \(\dot{x}\)
- Refinement:
  - same 2 points with derivatives

---

## Method

Construct a cubic interpolant matching both values and derivatives at interval endpoints.

This gives a smooth approximation for both:

\[
x(t)
\]

and

\[
\dot{x}(t)
\]

Threshold crossings are found by solving:

\[
x(t)=x^*
\]

Extrema are found by solving:

\[
\dot{x}(t)=0
\]

The extremum condition reduces to solving a quadratic equation.

---

## Pros

- Excellent accuracy/cost tradeoff
- Smooth interpolation
- Uses derivative information already available from ODE evaluation
- Accurate threshold and extremum refinement
- Standard approach in many ODE solvers

---

## Cons

- Slightly more implementation complexity
- Cubic root solve may require safeguarded iteration
- Extrema remain numerically ill-conditioned when very flat

---

# 4. Higher-order polynomial interpolation

Examples:
- quintic interpolation
- collocation polynomials
- solver-native dense output

---

## Buffer requirements

Typically:
- 3+ points
- often derivatives as well

---

## Pros

- Very high interpolation accuracy
- Can match high-order solvers

---

## Cons

- More complicated
- Greater oscillation risk
- Usually unnecessary for event localization
- Diminishing practical returns

---

# Event detection vs refinement

## Threshold crossing

### Detection

Usually based on a sign change:

\[
(x_0-x^*)(x_1-x^*) < 0
\]

Crossings are generally well-conditioned if:

\[
|\dot{x}| \not\approx 0
\]

Steep crossings localize accurately.

---

## Extremum

### Detection

Typically detected from:

\[
\dot{x}_0 \dot{x}_1 < 0
\]

or from curvature in 3-point state data.

### Conditioning

Extrema are harder because:

\[
\dot{x}=0
\]

Small interpolation errors can strongly perturb the estimated extremum time.

Broad, flat extrema are intrinsically difficult.

---

# Practical recommendation

For smooth ODE solutions:

## Recommended buffer

Store:

- \(t_n\)
- \(x_n\)
- \(\dot{x}_n\)

for the two most recent accepted solver steps.

---

## Recommended method

Use cubic Hermite interpolation for:

- threshold crossing refinement
- extremum refinement
- dense output

This is the best overall balance of:

- simplicity
- robustness
- smoothness
- computational efficiency
- accuracy

---

# References

## Cubic Hermite interpolation
- https://en.wikipedia.org/wiki/Cubic_Hermite_spline

## Hermite interpolation
- https://en.wikipedia.org/wiki/Hermite_interpolation

## Dense output for ODE solvers
- https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method

## Root finding
- https://en.wikipedia.org/wiki/Brent%27s_method

## Polynomial interpolation
- https://en.wikipedia.org/wiki/Polynomial_interpolation
