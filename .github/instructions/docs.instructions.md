---
applyTo: "docs/**"
description: Use when editing user docs under docs/. Keep docs workflow-shaped, current-package focused, and aligned with the docs IA plan.
---

# Docs Rules

- Keep docs centered on the current Python package, supported workflows, and reproducible usage.
- Prefer task-oriented guides over broad architectural narration.
- Keep `docs/index.md` as docs-site orientation, not a second copy of the package landing page.
- Treat `docs/examples.md` as a curated appendix, not the main place where concepts are introduced.
- Keep installation, getting-started, guides, examples, and reference roles distinct.
- Avoid development-history narration, removed internal layers, or migration archaeology in user-facing pages.
- If the docs information architecture, page ownership, or wording boundary is changing, consult `.design/reference/docs_layout_plan.md` and `.design/reference/public_surfaces_plan.md`.
- If the task coordinates wording across `README.md`, `docs/`, or `paper/paper.md`, use the `public-surface-sync` skill.
- Tie important numerical, continuation, and performance claims to runnable examples, docs pages, or maintained tests.