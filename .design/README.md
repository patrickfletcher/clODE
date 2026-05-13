# .design Context Guide

## Default Read Path

Stop as soon as you have enough context for the task.

1. `.design/package_state.md` for the live package map, code ownership, and current constraints.
2. `.design/next_pr.md` for the current narrow implementation target.
3. `.design/ideas.md` only if you need dependencies, adjacent work, or follow-on ideas.
4. `.design/development_roadmap.md` only if you need broader sequencing or architectural rationale.
5. `.design/reference/README.md` only for task-specific deep dives.
6. `.design/archived/README.md` only for historical rationale or bug archaeology.
7. `.design/tmp/README.md` only if a current PR explicitly points to scratch work there.

## What Is Authoritative

- The root `.design/*.md` files are the live planning surface and the default source of truth for current maintainer context.
- `.design/reference/` holds focused deep dives and audit notes. Read them only when the task needs their topic.
- `.design/archived/` is historical context only. It never overrides the live root docs.
- `.design/tmp/` is scratch space. Ignore it by default.
- Repo memory snapshots can be useful for quick context, but source-controlled `.design/` docs are the authoritative repo-local record.
- New implementation work should target canonical modules first, not compatibility barrels, unless the task is explicitly about compatibility cleanup.

## Minimum Read Sets

- Narrow implementation PR: `package_state.md` -> `next_pr.md`.
- Broader refactor or module move: `package_state.md` -> `development_roadmap.md` -> `ideas.md`.
- Tests or diagnostics: `package_state.md` -> `next_pr.md` -> `reference/README.md` -> the one matching reference note.
- Docs, paper, or repo-surface work: `reference/README.md` -> `project_principles.md`, `docs_layout_plan.md`, `public_surfaces_plan.md`, or `joss_audit.md`.
- Scope, positioning, or cross-cutting design decisions: `package_state.md` -> `development_roadmap.md` -> `reference/README.md` -> `project_principles.md`.
- Runtime/OpenCL internals: `package_state.md` -> `reference/README.md` -> `pyopencl_leverage_audit.md` or `continuation_timebase_note.md`.
- Historical investigation: current docs first, then `archived/README.md`.

## Live Docs At The Root

- `package_state.md`: factual repo and package map.
- `next_pr.md`: current active implementation target and acceptance criteria.
- `ideas.md`: terse backlog and dependency board.
- `development_roadmap.md`: longer rationale and priority ordering.
- `MAINTENANCE.md`: rules for adding, updating, moving, archiving, and pruning design docs.

## Editing Rules

- Update an existing live doc before creating a new one.
- Keep the root small. New deep dives should usually live under `.design/reference/`.
- Promote enduring facts into `package_state.md`; keep `next_pr.md` narrow and current.
- Move completed or stale notes into `.design/archived/` rather than leaving them at the root.
- Keep scratch reproductions in `.design/tmp/`, then delete, promote, or archive them before the PR is done.
- When the package layout, active target, or the meaning of a reference note changes, update the relevant `.design` doc in the same PR.
- Read `.design/MAINTENANCE.md` before reshaping this tree or adding a new design note.
