# Semantic Layout Background

Purpose: preserve the broader peer-comparison and package-layout option analysis removed from the slim live semantic-layout note.
Read when: you need contributor-facing layout rationale beyond the current live guidance or want the historical comparison set behind the concept-first recommendation.
Update when: a signpost to the current live note needs correction.

Historical only. The current live summary is `.design/reference/semantic_layout_audit.md`.

## Archived peer-library comparison set

This comparison set was useful for checking recurring concept boundaries across successful solver libraries, but it no longer needs to sit in the live note.

### `martenlienen/torchode`

Useful pattern that informed clODE:

- separate problem or term definitions from numerical methods, controller policy, solve wrappers, and returned solution objects

### `patrick-kidger/diffrax`

Useful pattern that informed clODE:

- keep the organization concept-first around terms, solver algorithms, controllers, events, output policy, solution objects, and randomness sources

Historical caution:

- Diffrax carries more abstraction depth than clODE needs today; the useful lesson was separation of concerns, not type-count parity

### `google-research/torchsde`

Useful pattern that informed clODE:

- make stochastic execution resources and extra solver state explicit instead of hiding them entirely in implementation detail

### `SciML/DifferentialEquations.jl`

Useful pattern that informed clODE:

- separate problem definition, algorithm choice, execution or integrator state, callbacks or events, and returned solution objects clearly

## Archived layout options from the broader audit

### Option A: incremental concept-first layout, keep kernels where they are

Shape:

- keep `clode.problem`, `clode.observers`, `clode.simulation`, `clode.runtime`, and `clode._opencl`
- add clearer Python-side semantic modules inside those packages, or a narrow `clode.steppers` package later
- keep `_opencl` focused on execution, runtime, build, and dispatch concerns
- keep `clode/kernels/` as the asset root for the time being

Why this was preferred:

- lowest churn
- preserves the current package map
- improves contributor navigation without forcing kernel relocation up front

### Option B: co-locate semantic Python definitions and OpenCL assets per concept

Why it was deferred:

- current source assembly assumes a single kernel root and centralized include tree
- moving kernel assets before semantic definitions stabilize risks turning conceptual cleanup into mostly path churn
- runtime and execution concerns would still need a separate home either way

### Option C: per-concept backend submodules

Why it was rejected for now:

- too much module churn for the current package size
- risks reintroducing framework-like indirection before the semantic models themselves are settled

## Historical near-term order from the broader audit

1. Introduce explicit IVP semantics with built-in batch shaping.
2. Revisit whether a dedicated ensemble type is warranted only after real batch-specific behavior appears.
3. Make per-work-item solver state explicit.
4. Introduce clearer stepper and observer definition objects.
5. Revisit kernel relocation only after those semantic models exist.
