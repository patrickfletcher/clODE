# Ideas Board Closeout 2026-05-29

Purpose: archive completed planning-board lines moved out of `.design/ideas.md` once their durable outcomes were reflected in the live design docs.
Read when: you need the original wording of completed `ideas.md` items or historical traceability for why open dependencies refer to already-landed work.
Update when: this archive description becomes inaccurate.

These lines were archived on 2026-05-29 after verifying that the live `.design` root or reference notes already carry the current durable outcome.

Current source of truth for live planning remains `package_state.md`, `next_pr.md`, `development_roadmap.md`, and the relevant files under `.design/reference/`.

## Observer And Feature Model

- [x] P1 Observer readout naming and legacy-family retirement: standardize feature-readout helper names and event-value labels across semantic observers so `summary`, `normalized_schmitt_trigger`, and `neighborhood_return` can fully replace the legacy `basic`/`basicall`, `threshold_2`, `neighborhood_2`, `threshold_1`, and `neighborhood_1` surfaces. depends: observer-specific parameter models/classes beyond the threshold-family seam. refs: `clode/observers/types.py`, `clode/observers/_definitions.py`, `clode/simulation/features.py`, `docs/feature_extraction.md`, `docs/observers.md`
