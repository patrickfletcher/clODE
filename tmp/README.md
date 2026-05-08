# tmp Planning Docs

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
