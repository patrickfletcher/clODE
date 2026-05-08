# tmp Planning Docs

## Start Here

Recommended read order for a fresh session:

1. `tmp/package_state.md`: current architecture, canonical module homes, active compatibility surfaces, and packaging/runtime constraints.
2. `tmp/next_pr.md`: the current narrow implementation target.
3. `tmp/ideas.md`: short living board for backlog items and sequencing.
4. `tmp/development_roadmap.md`: longer rationale and prioritization.
5. the focused audit notes only when the task needs them (`testing_audit.md`, `backend_strategy_audit.md`, `pyopencl_leverage_audit.md`, `module_layout_plan.md`, `opencl_boundary_audit.md`, `joss_audit.md`).

Repo memory should stay aligned with the same state:

- `/memories/repo/current-state.md`: concise live snapshot.
- `/memories/repo/build-notes.md`: build, packaging, and local runtime validation facts.

## What Is Authoritative

- `tmp/package_state.md` is the live package and repo map.
- `tmp/next_pr.md` is the live execution target.
- `tmp/ideas.md` is the short living backlog.
- `tmp/archived/` is historical context only. Do not treat it as the source of truth for the current package.
- New implementation work should target canonical modules first, not compatibility barrels, unless the task is explicitly about compatibility cleanup.

- `ideas.md`: short living planning board. Keep items terse, grouped by domain, and mark them done with checkboxes.
- `package_state.md`: factual map of what code lives where and what it currently does.
- `development_roadmap.md`: longer rationale and prioritization note.
- `module_layout_plan.md`: proposed role-based package layout after the PyOpenCL-first decision.
- `next_pr.md`: the current narrow implementation target.
- `joss_audit.md`: JOSS fit, positioning, and readiness gaps.
- `backend_strategy_audit.md`: whether another backend is realistic and what that implies for `_backends/`.
- `pyopencl_leverage_audit.md`: PyOpenCL features worth leveraging more aggressively.
- `testing_audit.md`: test taxonomy and kernel-component test strategy.
- `archived/`: historical migration notes, bug archaeology, and older design docs.

## Editing Rules

- Add new work to `ideas.md` first.
- Keep one-line items in `ideas.md` when possible; move detail into the roadmap, next-PR note, or a focused deep-dive doc only when needed.
- Use `depends:` for prerequisites and `blocks:` for work that should wait on the item.
- Prefer updating existing items over creating near-duplicates.
- If package structure, canonical import paths, or build/runtime assumptions change, update `tmp/package_state.md`, this index, and the relevant repo memory snapshot in the same pass.
