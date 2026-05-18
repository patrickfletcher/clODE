# Reference Notes

Purpose: task-specific deep dives and audits.
Read when: the root `.design` docs stop being enough and you need a focused note on testing, docs, runtime behavior, publication scope, or a specific design question.
Update when: the meaning of a topic changes or a note is promoted, split, or archived.

Read the root `.design` docs first. This directory is not default reading.

## Topics

- `project_principles.md`: settled project principles and open strategic design questions that should guide broad refactors, public wording, and scope decisions.
- `semantic_layout_audit.md`: current package-layout guidance for IVP-first batch semantics, simulator orchestration, solver state, stepper definitions, and observer definitions.
- `solver_state_implementation_plan.md`: code-facing implementation map for the active solver-state ownership PR, including file targets, phase ordering, and focused test coverage.
- `ivp_api_test_plan.md`: closeout note for the landed first-pass `InitialValueProblem` API shape, staged delivery order, Python-backed callability, and the slim test matrix.
- `continuation_timebase_note.md`: current solver-owned time-base semantics, attained-`tf` continuation guidance, and the remaining solver-state questions.
- `testing_audit.md`: test taxonomy and the kernel-component testing strategy.
- `pyopencl_leverage_audit.md`: PyOpenCL features and runtime helpers worth using more aggressively.
- `docs_layout_plan.md`: docs information architecture, example-execution policy, and API reference cleanup path.
- `public_surfaces_plan.md`: README, docs-home, package-index, and paper ownership and wording guardrails.
- `joss_audit.md`: publication and JOSS-readiness notes.

## Rules

- Do not treat these notes as overriding `package_state.md`, `next_pr.md`, `ideas.md`, or `development_roadmap.md`.
- If a reference note becomes active delivery scope, summarize the live decision back into the root docs.
- Archive or split a reference note once it stops being a useful current deep dive.
