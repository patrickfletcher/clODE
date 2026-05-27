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

- The next architectural gap is no longer the solver-owned failure-policy path. It is that integration settings, trajectory output policy, observer runtime settings, event-output capacity, persistent observer state, and fetched outputs are still represented across compatibility bundles, simulator subclasses, caches, and `_opencl` helpers rather than one explicit owner model.
- That ambiguity matters more than observer register pressure because it makes the code harder to reason about, blocks a cleaner observer-authoring story, and weakens the package's distinctive JOSS narrative around on-device observers and outputs.
- The package's numerical and continuation claims still need tighter exact-solution, convergence, and public evidence coverage, but that is now the second grouped follow-on rather than the immediate target.
- The heavier feature kernels still need an observer-state and register-pressure audit, but that is better treated as a later optimization pass after the semantic owner model is explicit.
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

### 1. Simulation state and output ownership model

Clarify which concepts are Python-owned semantic models and which are `_opencl` execution details: integration settings, trajectory output policy, observer runtime settings, event-output policy, solver state, persistent observer state, and fetched outputs.

Why first:

- it gives the codebase one legible answer to what solver state, observer state, and outputs mean on the Python side versus the OpenCL side
- it makes later observer authoring/composition, output ergonomics, and public narrative work easier to explain and safer to extend

### 2. Observer concept audit and authoring follow-through

Once the owner split is explicit, first audit what clODE observers currently mean, where that model overlaps with standard event-function semantics, and which parts of feature selection or memory-footprint control need separate declaration on the Python side. Then shape a narrower path for adding more built-in observers and, later, composable or user-authored observers from Python-side definitions instead of a hand-wired kernel catalog.

That strongest first proof target is now landed in the summary-observer family: `basic`, `basicall`, and `summary` now share one declaration family with explicit selection policy. The threshold-family follow-on proofs are now landed too: `Observer.threshold_crossing` proved that a useful absolute-threshold event family can keep threshold value and crossing direction in runtime settings while exposing only a lean timestamp-plus-count readout, `Observer.normalized_threshold_crossing` proved the same lean readout can coexist with warmup-derived normalized threshold parameterization, and the current Schmitt surface now spans both absolute and warmup-derived semantics through `Observer.schmitt_trigger` and `Observer.normalized_schmitt_trigger`, while `Observer.threshold_2` remains the retained legacy full normalized Schmitt path. That same semantic split is now landed for the adjacent trigger families too: `Observer.local_extremum` and `Observer.neighborhood_return` are the lean semantic observers, while `Observer.local_max` and `Observer.neighbourhood_2` remain retained heavier legacy workflows. The shared-config seam is now landed where it fits: the semantic observer names plus `ThresholdCrossingConfig`, `SchmittTriggerConfig`, `LocalExtremumConfig`, and `NeighborhoodReturnConfig` are preferred, while the broad parameter bundle and retained heavy observers remain compatible. The next useful move is to decide whether recurring oscillation-oriented controls, readouts, and state-size decisions deserve a shared bundle seam instead of further family-local duplication.

That same audit still suggests that `nhood1` should not guide the next abstraction unless a clearer deterministic workflow emerges. The landed `local_extremum` and `neighborhood_return` families answer the lean trigger-semantics question, while the more durable design question now is how much of the heavier oscillation-oriented readout bundle from `local_max`, `threshold_2`, and `nhood2` deserves a shared declaration seam at all. The retained heavy observers remain useful as evidence for that question, while slope-gated Schmitt behavior should still be framed as a narrower noisy-trace tool rather than the generic threshold model.

Why second:

- observers are a first-class clODE concept and one of the clearest ways to differentiate the package in both engineering and publication terms
- it avoided mixing higher-risk authoring-surface decisions into the lower-level owner-split PR, and the landed trigger-semantics slices kept the next event-observer step evidence-backed instead of speculative

### 3. Numerical evidence and publication follow-through

Keep the proof and publication story moving after the owner model is clearer: add a few more representative exact-solution problems, keep the release-gating evidence centered on `test/core_numerics/`, and build the benchmark or comparison material the paper will need.

Why third:

- it strengthens the JOSS case without letting public-surface work outrun the internal execution and observer model
- it keeps docs, benchmarks, and the paper tied to a clearer and more stable package story

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
