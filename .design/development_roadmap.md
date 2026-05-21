# clODE Development Roadmap

Purpose: medium-lived rationale and priority ordering for follow-on work.
Read when: you need tradeoffs, architectural sequencing, or the why behind a workstream.
Update when: priorities shift materially or the architectural rationale changes.

Use `.design/ideas.md` as the short living board. This file is the longer rationale and prioritization note.

## Current Read On The Repo

### Strengths

- The numerical core is still strong: one work-item per trajectory, compile-time specialization, and a clean separation between stepper kernels and observer kernels.
- The PyOpenCL migration put runtime ownership, source assembly, build keys, and program caching in Python where they are easier to reason about and test.
- `FeatureSimulator` remains a distinctive strength because it can extract useful oscillation and event statistics without forcing full trajectory storage.
- Observer definitions, resolved specs, stepper definitions, solver-state boundaries, and cache invalidation are now coherent enough that future work can target measured gaps rather than basic ownership cleanup.
- The repo already has a credible correctness backbone: exact-solution numerics, continuation regressions, and a small but useful `test/kernel_components/` layer.

### Friction That Still Matters

- The package's numerical and continuation claims are now ahead of the proof layer. The main near-term gap is tighter exact-solution, convergence, and public evidence coverage.
- The heavier feature kernels still need an observer-state and register-pressure audit before deeper redesign decisions are justified.
- Large-ensemble ergonomics are still missing the next obvious layer: IVP-side batch-generation helpers and device-capacity batching.
- Diverged-time continuation still has a real representational limit because the live shared-window model cannot encode exact continuation after different attained `tf` values.
- Public compatibility bundles still mix concerns: `SolverParams` spans integration and output policy, and `FeatureSimulator` still exposes broad legacy `observer_*` inputs.
- The runtime is intentionally single-device only, and there is still no good path for stiff or implicit methods.

### What Should Not Be Reopened Casually

- Keep the current semantic package split: `problem`, `observers`, `simulation`, `runtime`, and `_opencl`.
- Keep the root compatibility barrels for now; they still carry import-compatibility and collaborator-orientation value.
- Treat `InitialValueProblem`, first-pass solver state, observer definitions, stepper definitions, and the output-policy split as landed foundations.
- Treat shared-final-time continuation via `advance_tspan_to_attained_final_time()` as the current exact-continuation convenience path, and keep broader diverged-time work deferred.
- Treat the current fixed-step, adaptive-time, and observer-helper adoption slice as landed; future helper rollout should be evidence-driven rather than speculative.

## Recommended Near-Term Sequence

Keep the next planning pass centered on one product story: fast large-ensemble workflows where users mainly want final states or on-device features, and where single-precision robustness is part of the value rather than a later cleanup.

### 1. Numerical validation and evidence bundle

Keep the next validation pass narrow and purposeful: a slim exact-solution solver-validation suite, alignment of precision-sensitive references with the live solver semantics, and public docs or examples that show where the adopted safeguards matter.

Why first:

- it strengthens the proof layer without turning into a broad test-suite rewrite
- it lets public docs and planning notes point to maintained evidence instead of ad hoc demos

### 2. Observer-state and register-pressure audit

After the proof layer is tighter, audit the heavier feature kernels with throughput in mind. Measure per-work-item observer state, retained event storage, and avoidable private scratch before deciding whether deeper observer-time redesign is worth the complexity.

Why second:

- it targets the main remaining performance risk in large feature sweeps
- it keeps redesign pressure evidence-driven instead of speculative

### 3. Large-ensemble execution ergonomics

After the validation bundle and observer-footprint audit are clearer, group the next ergonomics work around IVP-side batch-generation helpers and device-capacity ensemble batching. Keep broader trajectory chunking, dense output, and output-surface expansion behind those more central workflows.

Why third:

- it aligns the next ergonomics work with the package's strongest workflow niche
- it avoids broad trajectory-surface work before final-state and feature-first workflows are fully hardened

### Later and lower leverage

Keep source-assembly reshaping, broader PyOpenCL helper adoption, compatibility-barrel cleanup, richer observer-surface expansion, diverged-time continuation follow-through, and implicit-method work behind the grouped passes above unless a concrete bug or benchmark result pulls one of them forward.

## Priority Guardrails

- Prefer proof, measurements, and concrete workflow wins over structural churn.
- Keep public API stability unless a change clearly pays for itself immediately.
- Treat README, docs, and paper work as current-package narrative, not implementation history.
- Use `.design/ideas.md` for one-line backlog items; keep this file for rationale and sequencing only.

## Longer-Term Items That Need Deliberate API Discussion

- true multi-device execution via a dedicated API rather than by extending the current selectors
- a public observer-definition or custom-observer interface
- explicit public solver-state or observer-state objects
- public implicit or linearly implicit steppers
- a cleaner public separation between integration parameters and output-policy parameters

## Packaging Note

- Omitting docs from the sdist is reasonable if artifact size ever becomes a real concern.
- Keeping tests in the sdist remains the better default because downstream maintainers benefit from having the verification surface available.
