# Ideas Board Closeout 2026-05-31

Purpose: archive completed planning-board lines moved out of `.design/ideas.md` once their durable outcomes were reflected in the live design docs.
Read when: you need the original wording of completed `ideas.md` items or historical traceability for why open dependencies refer to already-landed work.
Update when: this archive description becomes inaccurate.

These lines were archived on 2026-05-31 after verifying that the live `.design` root or reference notes already carry the current durable outcome.

Current source of truth for live planning remains `package_state.md`, `next_pr.md`, `development_roadmap.md`, and the relevant files under `.design/reference/`.

## OpenCL Program And Runtime Model

- [x] P1 Single-source-of-truth accepted-step `K`-sample solution-buffer or solver-state kernel abstraction plus shared update helpers for steppers and observers after the explicit solver-state model lands. depends: explicit solver-state model and settled observer solution-buffer audit. refs: `clode/kernels/observers.cl`, `clode/kernels/transient.cl`, `clode/kernels/features.cl`, `.design/reference/semantic_layout_audit.md`, `.design/reference/observer_solution_buffer_audit.md`, `.design/reference/observer_concept_audit.md`

  Outcome: `package_state.md` confirms `advanceAcceptedStepHistory2`, `advanceAcceptedStepHistory3`, and their `ByVariable` variants are live in `observers.cl` and consumed by all lean observer families. The `development_roadmap.md` records the solution-buffer audit as complete and the K-sample concept as confirmed. Remaining adoption by heavier observer families is carried forward under the "Oscillation-oriented observer bundles" and "Extremum-family follow-through" backlog items.

## Observer And Feature Model

- [x] P2 Retire or quarantine `nhood1` unless a deterministic workflow justifies keeping it as more than historical experimentation. depends: deterministic static-trigger family direction settled. refs: `clode/kernels/observers/observer_neighborhood_1.clh`, `.design/reference/observer_concept_audit.md`

  Outcome: this line was removed as stale. The current tree no longer contains `observer_neighborhood_1.clh`, and the public observer surface is now centered on `summary`, `threshold_crossing`, `normalized_threshold_crossing`, `schmitt_trigger`, `normalized_schmitt_trigger`, `local_max`, and `normalized_neighborhood_return`. Any future neighborhood-return follow-through is tracked under the normalized family items already on `ideas.md`.
