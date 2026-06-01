# clODE Development Roadmap

Purpose: medium-lived rationale and priority ordering for follow-on work.
Read when: you need tradeoffs, architectural sequencing, or the why behind a workstream.
Update when: priorities shift materially or the architectural rationale changes.

Use `.design/ideas.md` as the short living board. This file is the longer rationale and prioritization note.

## Current Read On The Repo

### Strengths

- The numerical core is strong: one work-item per trajectory, compile-time specialization, and a clean separation between stepper kernels and observer kernels.
- The PyOpenCL migration put runtime ownership, source assembly, build keys, and program caching in Python where they are easier to reason about and test.
- `FeatureSimulator` remains a distinctive strength because it can extract useful oscillation and event statistics without forcing full trajectory storage.
- Observer definitions, resolved specs, stepper definitions, solver-state boundaries, and cache invalidation are now coherent enough that future work can target measured gaps rather than basic ownership cleanup.
- The repo already has a credible correctness backbone: exact-solution numerics, continuation regressions, and a small but useful `test/kernel_components/` layer.

### Friction That Still Matters

- The owner split is much clearer than it was before the recent state/output pass, but observer-side modeling is still uneven: compatibility bundles remain at the public surface, the heavier observer families still carry family-local accepted-step buffers and private state/output conventions, and the next reuse boundary is not settled yet.
- That remaining ambiguity matters more than observer register pressure by itself because it blocks a cleaner observer-authoring and packaging story and keeps the next observer-bundle decision harder to reason about than it should be.
- The package's numerical and continuation claims still need tighter exact-solution, convergence, and public evidence coverage, but that is now the second grouped follow-on rather than the immediate target.
- The heavier feature kernels still need an observer-state and register-pressure audit, but that is better treated as a later optimization pass after the solution-buffer and bundle boundary is explicit.
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

Keep the next planning pass centered on one product story: clODE's differentiator is not just solving ODEs on OpenCL, but supporting first-class on-device observers and output workflows from Python with a model that remains inspectable and extensible.

### 1. Oscillation-oriented observer bundles

The solution-buffer audit is complete: a shared `K`-sample solution-buffer concept was confirmed as a prerequisite, K=2/K=3 shared accepted-step history update helpers are now live and adopted by all canonical event-observer families, and family event semantics remain local. The remaining near-term observer question is whether recurring event-trigger controls (`min_amp`, `max_event_count`) plus family-specific oscillation readouts should be factored through a shared seam or kept family-local.

Why next:

- observers are a first-class clODE concept and one of the clearest ways to differentiate the package in engineering and publication terms
- `max_event_count` is already a general event-loop limiter across the current event-triggering families, so documenting it as a canonical control rather than as legacy baggage will clarify both user docs and future observer authoring
- `min_amp` already has a subtle semantic difference between one-pass and warmup-derived families; naming it as an event gate on `event_var` variation will clarify that seam before any broader family adoption decision
- the current event families already share some combination of event-loop controls and oscillation-oriented readouts, so a clear shared seam would improve both maintainability and observer UX consistency without overstating derivative-gate fields such as `eps_dx`

### 2. Large-ensemble ergonomics: batch-generation helpers and device-capacity batching

Two open P1 items form a natural second pass once the observer-bundle direction is settled enough that the batch layer does not have to accommodate large structural observer changes: IVP-side helpers for grid, random, and quasi-random ensemble construction (currently buried in `Simulator`), and device-capacity batching for ensembles that exceed buffer limits.

Why second:

- large-ensemble use cases are one of clODE's stated strengths and both items are blocked on a clearer IVP ownership model that is already mostly in place
- the batch-generation helpers and the device-capacity batching share the same IVP-owned shape and output-policy prerequisites, so tackling them together avoids revisiting the same boundary twice
- landing this pass would make the quick-start story and example suite significantly stronger before the publication push

refs: `.design/ideas.md` → "Batch-generation helpers" and "Ensemble batching for device-capacity limits", `.design/reference/chunked_execution_audit.md`, `.design/reference/semantic_layout_audit.md`

### 3. Numerical evidence and publication follow-through

Keep the proof and publication story moving once the observer-buffer direction is clearer: add a few more representative exact-solution problems, keep the release-gating evidence centered on `test/core_numerics/`, and build the benchmark or comparison material the paper will need.

Why third:

- it strengthens the JOSS case without letting public-surface work outrun the internal execution and observer model
- it keeps docs, benchmarks, and the paper tied to a clearer and more stable package story
- the release gate is now explicit and scriptable through `tools/run_test_bundle.py` (`release` alias), so numerical/publication follow-through can target concrete evidence surfaces instead of ad-hoc test runs

### 4. Later and lower leverage

Keep observer-state/register-pressure cleanup, source-assembly reshaping, broader PyOpenCL helper adoption, compatibility-barrel cleanup, diverged-time continuation follow-through, and implicit-method work behind the grouped passes above unless a concrete bug or benchmark result pulls one of them forward.

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
