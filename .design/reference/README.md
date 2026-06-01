# Reference Notes

Purpose: task-specific deep dives and audits.
Read when: the root `.design` docs stop being enough and you need a focused note on testing, docs, runtime behavior, publication scope, or a specific design question.
Update when: the meaning of a topic changes or a note is promoted, split, or archived.

Read the root `.design` docs first. This directory is not default reading.

Long notes should be read in layers: start with `## Bottom line` and any `## Fast path` section before loading deeper sections.

## Routing By Topic

### Current design audits and boundary notes

- `project_principles.md`: settled project principles and open cross-cutting design questions.
- `compatibility_boundary_audit.md`: current canonical versus compatibility-only observer surfaces and the follow-on cleanup questions they create.
- `observer_concept_audit.md`: current observer semantics, current state/bundle question, and implementation-path guidance for observer authoring and memory-footprint control.
- `observer_readout_audit.md`: current observer readout contract, naming guardrails, and open seam decisions.
- `observer_readout_phases.md`: short readout-rollout status board and pointer to archived rollout history.
- `observer_event_feature_variable_contract.md`: family-by-family `e_var_ix`/`f_var_ix` routing contract.
- `observer_solution_buffer_audit.md`: accepted-step `K`-sample buffer audit, current prerequisite decision, and smallest follow-on proof target.
- `semantic_layout_audit.md`: stable package-layout and ownership guidance. Use for durable placement principles, not the active PR sequence.
- `solver_state_implementation_plan.md`: current solver-state boundary and the durable guardrails it leaves in place.
- `continuation_timebase_note.md`: solver-owned time-base semantics and attained-`tf` continuation guidance.
- `chunked_execution_audit.md`: future chunking and ensemble-batching direction plus the current guardrails that later chunking work should preserve.
- `pyopencl_leverage_audit.md`: PyOpenCL helpers and runtime/build opportunities.

### Stable guardrails for numerics and testing

- `single_precision_numerics_note.md`: verified float32 failure modes, current demos, and mitigation guardrails.
- `single_precision_ode_solver_guide.md`: compact float32 arithmetic, scaling, conditioning, and mixed-precision guidance for sensitive ODE and kernel work. Read when editing `clode/kernels/*` or other float32-sensitive solver math.
- `ode_event_interpolation_note.md`: interpolation and event-refinement tradeoffs for threshold and extremum outputs; use when choosing event-time or event-state refinement methods.
- `testing_audit.md`: test taxonomy, evidence strategy, and component-test direction.

### Public-surface plans

- `docs_layout_plan.md`: docs-site IA, example-execution policy, and API reference cleanup.
- `public_surfaces_plan.md`: README, docs-home, paper, contributor docs, and repo-metadata ownership.
- `joss_audit.md`: publication-readiness assessment and evidence gaps.

## Rules

- Prefer one note or one topical cluster over a full-directory sweep.
- Keep long notes skimmable: add a first-screen `## Fast path` or equivalent stop-guidance block before letting a note grow into a default full read.
- Keep note type obvious: stable guardrails should not double as active roadmap docs, and current design audits should archive rollout history once it stops guiding active work.
- If a note exceeds roughly 200 lines, treat that as a prompt to split or archive landed-history sections unless the entire note is still active decision support.
- Do not treat these notes as overriding `package_state.md`, `next_pr.md`, `ideas.md`, or `development_roadmap.md`.
- If a reference note becomes active delivery scope, summarize the live decision back into the root docs.
- Archive or split a reference note once it stops being a useful current deep dive.
