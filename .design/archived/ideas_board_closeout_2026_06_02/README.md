# Ideas Board Closeout 2026-06-02

Purpose: archive completed planning-board and inbox lines moved out of `.design/ideas.md` once their durable outcomes were reflected in the live design docs.
Read when: you need the original wording of completed observer-event-geometry or observer-time-window ideas, or historical traceability for why later backlog lines refer to those already-landed results.
Update when: this archive description becomes inaccurate.

These lines were archived on 2026-06-02 after verifying that the live `.design` root and reference docs already carry the current durable outcome.

Current source of truth for live planning remains `package_state.md`, `next_pr.md`, `development_roadmap.md`, and the relevant files under `.design/reference/`.

## Performance, Randomness, And Scaling

- [x] P2 Leaner observer time-window architecture for feature-rich event detectors: evaluate elapsed-only or local-`dt` sample geometry instead of carrying both absolute and elapsed buffered times everywhere, without reintroducing large-origin subtraction. depends: observer-state and register-pressure audit. refs: `clode/kernels/observers/observer_local_maximum.clh`, `clode/kernels/observers/observer_normalized_schmitt_trigger.clh`, `.design/reference/single_precision_numerics_note.md`

  Outcome: the semantic event observers now retain a window origin plus solver-relative elapsed history instead of mirrored absolute-time buffers, and the live storage contract is summarized in `.design/package_state.md` and `.design/reference/observer_solution_buffer_audit.md`.

## Inbox

- [x] Implement event state readouts alongside event time stamps to support full event geometry capture for all observers

  Outcome: the canonical event observers now expose retained event timestamps plus state/auxiliary geometry, and the live readout contract is summarized in `.design/package_state.md`, `.design/reference/observer_readout_audit.md`, and the user docs under `docs/feature_extraction.md` and `docs/observers.md`.