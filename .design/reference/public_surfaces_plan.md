# Public Surfaces Plan

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
- Keep history, migration rationale, and internal boundaries inside `.design/` and other maintainer-only notes.

## Recommended Surface Ownership

| Surface | Primary audience | Purpose | Should contain | Should avoid |
| --- | --- | --- | --- | --- |
| `README.md` | GitHub visitors and PyPI readers | Shared landing page and package index long description | one-paragraph positioning, install snippet, one quick-start example, capability summary, docs links, repo layout, license | redirects-only copy, deep API detail, migration history, long benchmark discussion, repo-relative links that break on package indexes |
| `docs/index.md` (target) | Docs readers | Entry page for the documentation site | task-oriented navigation, short package overview, where to start, links to install/tutorials/API/examples | full repo pitch, repeated badges, long scholarly framing |
| `docs/install.md` | Users installing clODE | Platform and runtime setup | Python support, OpenCL runtime expectations, verification steps, source install notes | broad package overview already covered elsewhere |
| `docs/getting_started.md` | New users | First successful run | one end-to-end example, core concepts, next docs links | benchmark claims, contributor guidance |
| tutorial/example pages under `docs/` | Users evaluating workflows | Task-oriented walkthroughs | runnable examples, plots, expected outputs, links back to API | duplicate installation or project overview copy |
| `docs/api_reference.md` | Existing users | Reference surface | API signatures and concise usage notes | introductory narrative |
| `paper/paper.md` | Reviewers and research readers | Scholarly framing | statement of need, workflow niche, comparison to alternatives, benchmarks, scientific use cases | install instructions, repo tour, maintainer notes |
| `CITATION.cff` | Researchers and citation tools | Machine-readable citation metadata | software citation, authors, DOI/version links when available | prose explanation |
| `CONTRIBUTING.md` | Contributors | Project contribution path | local setup, test/docs commands, issue/PR expectations, release-adjacent notes | user onboarding or solver overview |
| GitHub repo settings | Casual visitors | Fast project summary outside markdown | short description, topics, homepage/docs URL, social preview image | long text duplicated from `README.md` |

## Recommended Target State

1. `README.md` is the canonical GitHub landing page and the package index long description.
2. `docs/index.md` becomes the MkDocs home page.
3. `docs/index.md` owns the docs-home job.
4. `pyproject.toml` points to `README.md`, and `README.md` uses absolute URLs so it renders correctly on package indexes.
5. `paper/paper.md` is refreshed as a research-facing artifact, not a substitute for user documentation.

## Concrete Checklist

### Phase 1: Repo Landing Page

- [x] Replace the root `README.md` redirect with a standalone landing page.
- [x] Keep the wording focused on current user-facing workflows and supported inputs.
- [x] Include direct links to installation, getting started, feature extraction, trajectory simulation, and the API reference.
- [ ] Add a short citation section once `CITATION.cff` exists.

Citation metadata remains intentionally deferred until the paper and software citation details are stable enough to publish.

### Phase 2: Docs Site Split

- [x] Create `docs/index.md` as the dedicated docs homepage.
- [x] Update `mkdocs.yml` so `Home` points to `docs/index.md` instead of `docs/README.md`.
- [x] Narrow `docs/index.md` to orientation and navigation rather than package-index prose.
- [x] Keep the docs home page separate from the package long description.
- [x] Collapse the GitHub/PyPI split back to `README.md` as the single source of truth for the package overview.
- [x] Add a clear tutorials/examples landing page in `docs/` that curates the best scripts from `examples/`.
- [x] Add a performance-notes page that documents current benchmark scripts and reproducibility context without overclaiming.
- [x] Remove the overlapping fast-and-slow page rather than keeping two docs pages for the same trajectory-oriented workflow.
- [ ] Keep iterating on the docs IA rather than treating the current nav as final while API and example priorities are still shifting.
- [ ] Keep the examples surface intentionally narrow and let the guides carry most of the explanatory burden.
- [ ] Promote the performance-notes page to headline benchmark claims once they are supported by stable reproducible runs.

### Phase 3: Project Metadata

- [ ] Add `CITATION.cff`.
- [x] Add `CONTRIBUTING.md`.
- [ ] Decide whether to add `CODE_OF_CONDUCT.md`.
- [x] Expand `[project.urls]` in `pyproject.toml` with Issues and, once ready, Changelog or Paper links.
- [ ] Fill out the GitHub About section, topics, and social preview image.

`CONTRIBUTING.md` now covers local setup, test bundles, docs builds, and the current public-surface ownership model. GitHub About metadata remains manual repo-settings work.

### Phase 4: Paper Refresh

- [x] Update the paper title and summary so they describe the current Python package directly.
- [x] Rework the statement of need around large ensemble simulation, online feature extraction, and model-ingestion workflows.
- [x] Refresh the alternatives section with present-day comparisons.
- [ ] Add benchmark and reproducibility material that can be rerun from the repo.
- [x] Align paper wording with the public docs so the project tells one consistent story.

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
