---
name: public-surface-sync
description: 'Align public-facing surfaces. Use when the user asks to refresh README.md, docs/, `paper/paper.md`, package-index copy, or public wording so the repository tells one current, evidence-backed story.'
argument-hint: '[optional: README, docs, paper, or a wording/topic area to focus on]'
user-invocable: true
---

# Public Surface Sync

Keep README, docs pages, the paper, and related public wording aligned without collapsing their distinct jobs.

## When To Use

- The user asks to align or refresh `README.md`, `docs/`, and `paper/paper.md`.
- The user wants public-facing wording audited against the current package state.
- The user wants the package narrative, docs-home copy, or paper framing updated together.
- A feature or workflow change means multiple public surfaces now need the same current-package story.

## Required Reads

Always read:

1. `.design/reference/public_surfaces_plan.md`

Then add only what the scope needs:

- docs pages or docs navigation: `.design/reference/docs_layout_plan.md`
- paper, publication claims, or citation-adjacent wording: `.design/reference/joss_audit.md`
- current package facts or capability boundaries: `.design/package_state.md`
- touched public files and only the smallest supporting docs, examples, or tests needed to back the claims you are editing

Use `design-doc-audit` instead when the task is about `.design` itself rather than public-facing copy.

## Default Bias

- Keep each public surface doing its own job.
- Describe the current Python package and supported workflows directly.
- Avoid development-history narration, removed internal layers, and implementation archaeology.
- Do not widen claims beyond the current evidence base.
- Prefer the smallest coordinated set of edits that keeps the public story consistent.

## Procedure

1. Identify which public surfaces are in scope and what claims or wording need to stay aligned.
2. Load the minimum evidence set needed to support those claims.
3. Update each touched surface according to its role: landing page, docs-home or guide, contributor note, or research-facing paper.
4. If the edits change the meaning of a public-surface policy note, update `.design/reference/public_surfaces_plan.md` or the smallest related note in the same pass.
5. Validate the touched markdown and config files.

## Success Criteria

- README, docs, and paper tell one current, evidence-backed story.
- `README.md` stays concise and package-index safe.
- The docs site stays task-oriented rather than duplicating the landing page.
- The paper stays research-facing rather than turning into user onboarding.