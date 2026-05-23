---
description: Use for all clODE work. Route to the smallest relevant repo-local design context first.
---

# Start Here

Use `.design/README.md` as the repo-local router.

- For trivial or localized fixes, direct follow-ups, or test reruns that do not change package semantics, public docs, or planning, read `.design/README.md` and stop there unless you need more context.
- For package, docs, testing-strategy, or planning changes, follow the minimum-path read order in `.design/README.md` and stop once you have enough context.
- Prefer file-scoped instructions and the smallest relevant `.design` note over broad repo sweeps.
- Repo memory can be a shortcut, but source-controlled `.design` docs win if they disagree.

# Skills

- Use the `design-ideas-inbox` skill when routing `.design/ideas.md` `## Inbox` items.
- Use the `design-doc-audit` skill when checking `.design` for stale priorities, blockers, duplication, or drift.

# Repo Rules

- New implementation work belongs in `clode.problem`, `clode.observers`, `clode.simulation`, `clode.runtime`, and `clode._opencl` unless the task is explicitly about compatibility behavior.
- Treat the flat root modules as compatibility barrels, not the default home for new implementation work.
- When touching OpenCL code under `clode/kernels/`, read `.design/reference/single_precision_ode_solver_guide.md` first, then add the smallest repo-specific numerics note needed for the task.
- When package layout, public behavior, the active implementation target, or the meaning of a reference note changes, update the smallest authoritative `.design` doc in the same PR.
- Public-facing docs should describe current behavior, supported workflows, and tradeoffs without development-history narration.

# python environment

Use the `clode` python environment for all clODE work:
```bash
~/envs/clode/bin/python
```