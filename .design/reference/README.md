# Reference Notes

Purpose: task-specific deep dives and audits.
Read when: the root `.design` docs stop being enough and you need a focused note on testing, docs, runtime behavior, publication scope, or a specific design question.
Update when: the meaning of a topic changes or a note is promoted, split, or archived.

Read the root `.design` docs first. This directory is not default reading.

## Routing By Topic

### Package architecture and execution

- `project_principles.md`: settled project principles and open cross-cutting design questions.
- `semantic_layout_audit.md`: package-layout guidance for IVP-first batch semantics, simulator orchestration, solver state, stepper definitions, and observer definitions.
- `solver_state_implementation_plan.md`: current solver-state boundary and what the landed first pass deferred.
- `continuation_timebase_note.md`: solver-owned time-base semantics and attained-`tf` continuation guidance.
- `chunked_execution_audit.md`: future chunking and ensemble-batching direction.
- `pyopencl_leverage_audit.md`: PyOpenCL helpers and runtime/build opportunities.

### Numerics and testing

- `single_precision_numerics_note.md`: verified float32 failure modes, current demos, and mitigation guardrails.
- `testing_audit.md`: test taxonomy, evidence strategy, and component-test direction.

### Public docs and publication

- `docs_layout_plan.md`: docs-site IA, example-execution policy, and API reference cleanup.
- `public_surfaces_plan.md`: README, docs-home, paper, contributor docs, and repo-metadata ownership.
- `joss_audit.md`: publication-readiness assessment and evidence gaps.

## Rules

- Prefer one note or one topical cluster over a full-directory sweep.
- Do not treat these notes as overriding `package_state.md`, `next_pr.md`, `ideas.md`, or `development_roadmap.md`.
- If a reference note becomes active delivery scope, summarize the live decision back into the root docs.
- Archive or split a reference note once it stops being a useful current deep dive.
