# clODE Project Principles

Purpose: canonical home for settled project principles and open cross-cutting design questions.
Read when: you are making a broad refactor, changing public wording, proposing new solver or observer capabilities, or deciding whether work belongs in clODE at all.
Update when: a cross-cutting principle changes, an open strategic question is settled, or the public narrative constraints materially shift.

This note is not a backlog and does not override the root `.design` docs.

- Use `.design/package_state.md` for current package facts and module ownership.
- Use `.design/next_pr.md` and `.design/ideas.md` for active work.
- Use this note for the durable "what kind of project clODE is" layer.

## Settled Principles

### Numerical And Scientific Priorities

- Preserve numerical behavior and correctness before simplifying internals or broadening scope.
- Treat the current OpenCL ensemble execution model as core to the package: one work-item advances one ODE instance, builds are specialized by precision, stepper, observer, and problem shape, and Python-side array layout remains Fortran-order unless there is an explicit design decision to change it.
- Keep feature extraction as a first-class workflow, not a side effect of trajectory simulation. Stateful observers and on-device feature or event computation are part of clODE's distinctive value.
- Treat robust single-precision mitigation as part of the product value for long-window and large-origin ensemble workflows, and widen mitigation scope only when the empirical evidence stays clear.
- Prefer exact-solution, contract, and parity-style regressions when changing numerical or semantic behavior.
- Favor explicit continuation and reproducibility semantics over convenience heuristics that hide the difference between requested time windows and attained final time.
- Prefer reusable kernel-side numerical helpers and component-level regressions for precision-sensitive observer and stepper logic over duplicated ad hoc arithmetic hidden inside many individual kernels.

### Workflow And Ecosystem Positioning

- Optimize for large ensembles, parameter sweeps, repeated stochastic realizations, and online feature extraction from Python.
- Center the public narrative on large-ensemble workflows where the desired output is a final state or an on-device feature or event summary; full trajectory output stays supported, but it is not the defining center of the package.
- Lean capability growth, docs, and examples toward dynamical-systems and bifurcation-analysis workflows when that strengthens clODE's large-ensemble niche without weakening the numerical story.
- Support multiple model-ingestion paths when they serve the same core workflow: typed Python RHS functions, OpenCL source files, and XPP models.
- Keep the package narrative workflow-shaped rather than trying to present clODE as the broadest general-purpose ODE ecosystem.
- Position clODE as a pragmatic niche between CPU-first general-purpose integrators, broader Julia solver ecosystems, XPPAUT-style model-authoring workflows, and ML-oriented differentiable ODE stacks.
- Keep performance claims concrete, reproducible, and tied to workload shape rather than vague accelerator rhetoric.

### Public API And Package Surface

- `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` are the durable user-facing entry points, with `clode` as the stable top-level import surface.
- Keep public APIs workflow-shaped around IVPs, simulators, observers, trajectories, and features rather than exposing source bundles, transfer buffers, or cache objects as first-class concepts.
- Runtime selection is explicitly single-device today. If multi-device execution is ever pursued, it should use a dedicated API rather than overloading the current selection arguments.
- New implementation work should land in `clode.problem`, `clode.observers`, `clode.simulation`, `clode.runtime`, and `clode._opencl`, not in the flat compatibility barrels.
- The flat compatibility barrels exist for import-path continuity, downstream stability, transition safety, and collaborator orientation while the semantic layout settles. They are not preferred homes for new code.
- Public-facing docs, examples, and the paper should describe the current Python package and supported workflows, not removed wrappers, migration history, or internal archaeology.
- Public API cleanup should happen through deliberate deprecation and documentation shifts rather than incidental breakage during internal refactors.

### Refactoring End State

- The end state of the current refactoring is not a brand-new public API. It is a stable internal semantic model underneath the existing workflow-shaped public surface, with thin compatibility layers kept close to the legacy user-facing entry points until deliberate deprecation work begins.
- Python should own the meaning of the core concepts: IVP/problem data, solver state, observer definitions, stepper definitions, output policy, and compile-time build specification. `_opencl` should consume those definitions to execute efficiently, not define them first.
- Each concern should have one clear semantic owner: IVP for next-solve problem data, solver state for execution progress and continuation facts, observer definitions for feature schema and warmup/event semantics, persistent observer state for cross-solve observer continuity, output policy for retained trajectory/event data, and build specification for compile-time specialization inputs.
- Program rebuilds, buffer reallocations, and runtime cache invalidation should follow explicit policy boundaries rather than broad mirrored-state resets.
- Chunking, batching, continuation helpers, and later public ergonomics should be orchestration layers over those stable owners, not new semantic owners themselves.
- Broader public config redesign should wait until observer-definition and stepper-definition work have landed, so the eventual public cleanup reflects stable internal boundaries rather than another moving target.

### Semantic Ownership And Internal Modeling

- Keep one semantic owner per concern. IVP-owned problem data, solver-owned execution state, persistent observer state, fetched outputs, and compile-time build specification should not share ownership.
- Separate compile-time kernel-specialization inputs from runtime state. Program rebuilds and cache keys should follow an explicit build specification rather than incidental cache invalidation or mirrored host fields.
- Treat device buffers and host mirrors as implementation details of those semantic owners. They should implement a semantic contract, not define solver or observer state.
- Prefer Python-owned semantic definitions and OpenCL-owned execution. Stepper, observer, and state meaning should be legible in Python first, with `_opencl` consuming those definitions efficiently.
- Prefer small, value-oriented, testable abstractions over backend-agnostic framework layers; do not generalize for hypothetical future backends before a second real backend exists.

### Implementation And Maintenance Guardrails

- Prefer Python-owned runtime, source-assembly, and metadata logic when it improves inspectability and maintenance without weakening numerical behavior.
- Introduce the smallest useful abstraction layer. Avoid growing framework-like indirection before the underlying solver and observer semantics are clear.
- Before broadening solver families or reorganizing kernel files, add component-level tests and shared helper layers where they reduce ambiguity in the current execution path.
- Solve correctness and state-model debt before expanding scope into broader solver families, larger API surfaces, or nominal multi-device features.
- Keep the root `.design/` surface small and authoritative; use reference notes for focused deep dives and cross-cutting rationale.

## Open Strategic Decisions

These are active design questions, not commitments.

### State Semantics And Public Ergonomics

- How explicit should solver state, observer state, RNG state, requested `t_span`, and attained final time become internally, and should any of that become public API?
- If clODE eventually gains a first-class continuation helper or API, should that wait until per-work-item time ownership is explicit enough that divergent-`tf` ensembles are not misrepresented?
- How much shared numerical-helper infrastructure should live in `clode/kernels/clODE_utilities.cl` and related includes before more ambitious observer or stepper expansion begins?

### Solver And Observer Architecture

- How far should integration policy be separated from output and storage policy (`max_store`, `nout`, event capacity) in the internal model and, later, in the public API?
- How explicit should the observer-definition model become, and should custom or composable observers eventually be user-facing?
- How far should clODE extend into stiff or implicit solvers and Jacobian generation without losing its current workflow focus?

### Scope And Platform Decisions

- Should the flat compatibility barrels eventually be deprecated, leaving `clode` as the primary stable top-level barrel?
- Which adjacent capabilities belong in clODE itself versus a sibling or helper package: broader interop, model-conversion tooling, benchmarking helpers, or analysis add-ons?

### Narrative And Evidence Thresholds

- What benchmark and comparison package is sufficient before making stronger public claims in the docs and paper?
- How much of this maintainer-facing principles layer should eventually be promoted into the public narrative versus kept as internal guidance?
