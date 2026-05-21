# Public Surfaces Plan

Purpose: define the distinct jobs and wording guardrails for README, the docs site, the package index page, the paper, contributor docs, and repo metadata.
Read when: editing README, the paper, contributor-facing docs, or cross-surface public wording and ownership.
Update when: public-surface ownership, wording constraints, or publication-adjacent repo metadata priorities change.

Use `docs_layout_plan.md` for docs-site IA and API-reference curation. This note is about cross-surface ownership and wording.

## Goal

- Give the repository, docs site, package index page, and paper distinct jobs.
- Keep public wording centered on the current Python package and supported user workflows.
- Avoid public references to removed internal layers, migration history, or implementation archaeology.
- Refine the docs information architecture incrementally while the public API and examples continue to evolve.

## Wording Guardrails

- Say `Python package`, `OpenCL-capable CPUs and GPUs`, `ensemble simulation`, `parameter sweeps`, `feature extraction`, and `trajectory simulation`.
- Describe supported inputs in user terms: Python RHS functions, OpenCL source files, and XPP models.
- Keep performance claims concrete and workload-shaped rather than absolute.
- Do not mention legacy wrappers, removed middle layers, or internal rewrites in public-facing copy.
- Do not narrate public features through development-history framing such as `old vs current`, `previously`, or `after the refactor`; describe the current package behavior and compare it to alternative algorithms only when that comparison helps users understand the current feature.
- Keep history, migration rationale, and internal boundaries inside `.design/` and other maintainer-only notes.

## Recommended Surface Ownership

| Surface | Primary audience | Purpose | Should contain | Should avoid |
| --- | --- | --- | --- | --- |
| `README.md` | GitHub visitors and PyPI readers | Shared landing page and package index long description | one-paragraph positioning, install snippet, one quick-start example, capability summary, docs links, repo layout, license | redirects-only copy, deep API detail, migration history, long benchmark discussion, repo-relative links that break on package indexes |
| `docs/index.md` | Docs readers | Entry page for the documentation site | task-oriented navigation, short package overview, where to start, links to install/tutorials/API/examples | full repo pitch, repeated badges, long scholarly framing |
| `docs/install.md` | Users installing clODE | Platform and runtime setup | Python support, OpenCL runtime expectations, verification steps, source install notes | broad package overview already covered elsewhere |
| `docs/getting_started.md` | New users | First successful run | one end-to-end example, core concepts, next docs links | benchmark claims, contributor guidance |
| tutorial/example pages under `docs/` | Users evaluating workflows | Task-oriented walkthroughs | runnable examples, plots, expected outputs, links back to API | duplicate installation or project overview copy |
| `docs/api_reference.md` | Existing users | Reference surface | API signatures and concise usage notes | introductory narrative |
| `paper/paper.md` | Reviewers and research readers | Scholarly framing | statement of need, workflow niche, comparison to alternatives, benchmarks, scientific use cases | install instructions, repo tour, maintainer notes |
| `CITATION.cff` | Researchers and citation tools | Machine-readable citation metadata | software citation, authors, DOI/version links when available | prose explanation |
| `CONTRIBUTING.md` | Contributors | Project contribution path | local setup, test/docs commands, issue/PR expectations, release-adjacent notes | user onboarding or solver overview |
| GitHub repo settings | Casual visitors | Fast project summary outside markdown | short description, topics, homepage/docs URL, social preview image | long text duplicated from `README.md` |

## Current Target State

1. `README.md` is the canonical GitHub landing page and the package index long description.
2. `docs/index.md` owns the docs-home job.
3. The docs site keeps a task-oriented IA and a narrow curated examples surface.
4. `paper/paper.md` remains a research-facing artifact rather than a substitute for user documentation.
5. `CONTRIBUTING.md` remains contributor-facing rather than part of the user docs narrative.

Most of the initial surface split is already landed. This note should track current policy and remaining open work, not a closeout checklist.

## Open Follow-Ups

- Add `CITATION.cff` once the paper and software citation details are stable enough to publish.
- Decide whether to add `CODE_OF_CONDUCT.md`.
- Fill out the GitHub About section, topics, and social preview image manually in repo settings.
- Keep iterating on the docs IA while the public API and example priorities continue to move.
- Keep the examples surface intentionally narrow and let the guides carry most of the explanatory burden.
- Add benchmark and reproducibility material that can be rerun from the repo before making stronger public performance claims.

## Example-source policy

- For short, deterministic flagship examples, it can make sense to keep a runnable script under `examples/` and include that file in the docs.
- The current lightweight path uses `pymdownx.snippets` plus a `py source run` fence so the docs can show Python source and render output below it.
- Do not force every docs snippet into this path. Keep inline snippets inline when the docs need a shorter teaching example than the repository script.
- For now, keep the examples surface selective: favor guide-backed examples and a small curated examples page over a large gallery.
- The broader docs information architecture and API reference cleanup recommendations live in `.design/reference/docs_layout_plan.md`.

## Deferred items

- `CITATION.cff`, `CODE_OF_CONDUCT.md`, and manual GitHub repo metadata are intentionally deferred for now.

## Notes On Repo Hygiene

- `site/` and `clode.egg-info/` are generated local outputs and should remain out of the maintained source narrative.
- Public pages should link to `docs/`, `examples/`, and `paper/` intentionally; generated output directories should stay ignored and unpublished in source control.
