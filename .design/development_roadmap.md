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

- The next architectural gap is not the low-level time-base itself; it is the remaining solver-owned per-work-item state follow-through for current-time and other unresolved time-base facts beyond the landed status, accepted-step-count, last-accepted-step-width, and collapsed-window no-progress path.
- The package's numerical and continuation claims still need tighter exact-solution, convergence, and public evidence coverage, but that is now the second grouped follow-on rather than the immediate target.
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

### 1. Solver-owned per-work-item state and failure reporting

The next grouped pass should extend the landed solver-owned status, accepted-step-count, last-accepted-step-width, and collapsed-window no-progress path to any remaining current-time values and any broader precision-loss reporting, without reintroducing solver-owned diagnostics through observer bookkeeping or observer public outputs.

Why first:

- the current wrapper boundary already preserves and surfaces accepted step width, accepted step counts, and failure status, so the remaining work is mainly about extending that ownership model to the unresolved time-base facts
- it restores single source of truth before any broader observer-memory or continuation redesign

### 2. Numerical validation and evidence bundle

After the solver-state boundary is clearer, keep the validation pass narrow and purposeful: a slim exact-solution solver-validation suite, alignment of precision-sensitive references with the live solver semantics, and public docs or examples that show where the adopted safeguards matter.

Why second:

- it strengthens the proof layer without mixing documentation work into the state-ownership refactor
- it lets public docs and planning notes point to maintained evidence instead of ad hoc demos

### 3. Observer-state and register-pressure audit

After the solver-state split is explicit and the proof layer is tighter, audit the heavier feature kernels with throughput in mind. Measure per-work-item observer state, retained event storage, and avoidable private scratch before deciding whether deeper observer-time redesign is worth the complexity.

Why third:

- it avoids auditing observer footprint while solver-owned diagnostics are still mixed into observer state
- it keeps redesign pressure evidence-driven instead of speculative

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
