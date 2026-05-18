# Docs Layout Plan

## Goals

- Keep the docs slim, task-oriented, and easy to scan.
- Prefer incremental improvements over a full docs rewrite while the public API and examples are still settling.
- Separate tutorials, how-to guides, examples, and API reference instead of mixing them on the same page.
- Use automatic example execution selectively for short, deterministic flagship examples only.
- Make the API reference readable by curating what is shown and by improving docstrings on the public types that actually matter.

## Current issues

- The docs nav is still page-oriented rather than task-oriented.
- The docs structure should keep evolving iteratively instead of being treated as locked before the API surface stabilizes.
- Several pages mix concepts, examples, and reference details.
- The current API reference page is noisy because it is driven by a broad mkdocstrings dump rather than a curated public reference.
- The observer-visualization examples exist as scripts, but their explanatory role should be clarified inside the features/observers guidance rather than left as a disconnected example cluster.

## Current iteration target

### Home

- `docs/index.md`: short orientation page with links into the main workflows.

### Getting started

- `docs/install.md`
- `docs/getting_started.md`
- `docs/querying_opencl.md` for runtime and device inspection
- simulator orientation should stay close to the getting-started material rather than becoming a large separate section too early

### Guides

- `docs/init_runtime.md`
- `docs/specifying_odes.md` for Python, OpenCL, and the overall model-definition path
- `docs/xpp_files.md` for XPP-specific ingestion details
- `docs/feature_extraction.md` as the main feature-and-observer guide
- `docs/trajectory_simulation.md`
- `docs/Ornstein-Uhlenbeck_process.md`
- runtime-specific guidance should stay narrow and practical rather than becoming a separate architecture narrative

### Features and observers

- Keep `docs/observers.md` as supporting concept material for the feature-extraction workflow.
- Treat `examples/visualize_events_localmax.py`, `examples/visualize_events_nhood2.py`, and `examples/visualize_events_threshold2.py` as the natural source material for explaining the main built-in observers with plots.
- Prefer folding that explanation into the features/observers guidance instead of building a separate top-level observer silo unless the content grows materially.

### Examples

- Keep `docs/examples.md` intentionally narrow for now.
- Lead with `examples/spike_counting.py` and add other example pages only when they illuminate a workflow that is not already covered well by a guide.
- Treat the examples surface as a curated appendix, not as the primary place where concepts are introduced.

### Reference

- `docs/api_reference.md`
- `docs/logging_levels.md`
- `docs/performance_notes.md`

### Out of docs nav

- `CONTRIBUTING.md` should stay contributor-facing rather than appearing in the user docs nav.
- The paper remains a separate research-facing artifact, not a docs section.

## Working nav direction

Do not treat this as locked. It is the current iteration target, not a promise that the final docs IA is settled.

- Home
- Getting Started
  - Installation
  - Getting Started
- Guides
  - Specifying ODEs
  - Runtime and device inspection
  - Choosing a simulator
  - Features and observers
  - Trajectory Simulation
  - Stochastic Simulation
- Examples
  - Spike Counting
  - Additional curated examples later, if needed
- Reference
  - API Reference
  - Logging and diagnostics
  - Performance Notes
- License

## Example execution policy

- Auto-run only short, deterministic examples that are central to the docs narrative.
- Prefer keeping heavyweight, exploratory, or device-sensitive scripts in `examples/` without running them during the docs build.
- For single-sourced examples, store the script in `examples/` and include it into the docs page only when the page genuinely benefits from showing the whole script.
- Keep short explanatory snippets inline when the full repository script would add noise.

## API reference recommendations

### Immediate fixes

- Use a curated API reference page instead of dumping private, inherited, and undocumented members.
- Group the reference by public concepts: simulators, outputs, configuration types, runtime query, and model-conversion helpers.
- Hide source listings by default in the rendered API reference.
- Render visible root headings for documented symbols and align mkdocstrings heading levels with the surrounding page sections so object titles and member titles use a predictable hierarchy.

### Docstring strategy

- Use one consistent public docstring style. Google-style docstrings are a good fit with mkdocstrings here.
- Start with the types that users actually inspect in the reference:
  - `Simulator`, `FeatureSimulator`, `TrajectorySimulator`
  - `ObserverOutput`, `TrajectoryOutput`
  - `InitialValueProblem`, `SolverParams`, `ObserverParams`
  - `DeviceInfo`, `PlatformInfo`, `query_opencl`, `print_opencl`
- Do not try to document every compatibility barrel or every exported builtin at once.

### Suggested next docstring targets

- Constructor and main method docstrings for the three simulator classes.
- Result-object method docstrings that explain returned shapes and structured-array behavior.
- Observer and solver parameter dataclasses so users can understand the knobs without reading source.

## Recommendation

- Reorganize the docs around tasks and workflows first, but do it incrementally.
- Keep the automatic example-execution path as a selective tool, not as a rule for every page.
- Keep the examples surface small until there is a clearer need for a broader gallery.
- Treat API reference cleanup as a mix of configuration cleanup and targeted docstring work, not as a pure formatting issue.
