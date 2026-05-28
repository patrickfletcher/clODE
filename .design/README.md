# .design Context Guide

Use this file as a router, not as background reading. Load the smallest relevant set and stop.

## Stop Rules

- If the task is a trivial or localized code or test fix and does not change package semantics, public docs, or planning, stop after this file and inspect the touched code.
- If the task changes package behavior, layout, public docs, testing evidence, or planning, keep reading only the matching minimum set below.
- Do not scan `.design/reference/`, `.design/archived/`, or `.design/tmp/` unless the task needs that topic.

## Minimum Read Sets

- Local bugfix or test repair: this file only, then code.
- Narrow implementation PR: `package_state.md` -> `next_pr.md`.
- Broader refactor or module move: `package_state.md` -> `development_roadmap.md` -> `ideas.md`.
- Tests, diagnostics, or numerical evidence: `package_state.md` -> `next_pr.md` -> `reference/README.md` -> the matching note.
- Docs, README, paper, or public-surface work: `reference/README.md` -> `project_principles.md`, `docs_layout_plan.md`, `public_surfaces_plan.md`, or `joss_audit.md`.
- `.design` maintenance or planning cleanup: `MAINTENANCE.md` -> `package_state.md` -> `next_pr.md` -> `ideas.md`.
- Historical investigation: current docs first, then `archived/README.md`.

## What Is Authoritative

- The root `.design/*.md` files are the live planning surface and the default source of truth for current maintainer context.
- `.design/reference/` holds focused deep dives and audit notes. Read them only when the task needs their topic.
- `.design/archived/` is historical context only. It never overrides the live root docs.
- `.design/tmp/` is scratch space. Ignore it by default.
- Repo memory can be a shortcut, but source-controlled `.design/` docs are the authoritative repo-local record.

## Root File Roles

- `package_state.md`: factual current package map, constraints, and contributor routing.
- `next_pr.md`: one active implementation target with acceptance criteria.
- `ideas.md`: terse open backlog and dependency board; completed lines belong in `.design/archived/` after the live docs reflect the landed result.
- `development_roadmap.md`: medium-lived rationale and sequencing.
- `MAINTENANCE.md`: rules for keeping this tree small and current.

## Editing Rules

- Update the smallest authoritative doc instead of creating a near-duplicate.
- Keep the root small. New deep dives should usually live under `.design/reference/`.
- Promote enduring facts into `package_state.md`; keep `next_pr.md` narrow and current.
- Move completed or stale notes into `.design/archived/` rather than leaving them at the root; `ideas.md` should not be used as completion history.
- Keep scratch reproductions in `.design/tmp/`, then delete, promote, or archive them before the PR is done.
- When the package layout, active target, or the meaning of a reference note changes, update the relevant `.design` doc in the same PR.
- Read `.design/MAINTENANCE.md` before reshaping this tree or adding a new design note.
