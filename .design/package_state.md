# clODE Package State

Purpose: factual package and repo map for the live codebase.
Read when: you need to know where code belongs, which public surfaces are live, or which constraints still apply.
Update when: canonical module homes, public compatibility surfaces, packaging rules, or runtime assumptions change.

## Snapshot

- clODE is now a pure-Python, PyOpenCL-only package.
- The public API is centered on `Simulator`, `TrajectorySimulator`, and `FeatureSimulator`, but the simulator classes are increasingly better understood as orchestration objects rather than the semantic owners of ensembles or solver state.
- Public simulators now construct the `_opencl` executors directly; the migration-era backend protocol/factory facade is gone.
- Split-window continuation correctness is now covered by live regressions for `basicall` features and seeded stochastic Euler on the current PyOpenCL path.
- Runtime logging now uses standard Python logging via `clode.configure_logging(...)` and `clode.get_logger(...)`; the old log-level compatibility API has been removed.
- Runtime-critical OpenCL source assets live under `clode/kernels/` and ship as package data.
- Historical migration and bug notes are archived under `.design/archived/`; they remain useful for rationale, but they are not the source of truth for the live package.
- A first-pass `InitialValueProblem` now owns default state, default parameters, basic batch shaping, remembered ensemble shape, and Python-backed SciPy-style RHS callability at the simulator boundary.
- Lower-level helper types such as `ProblemInfo` and `RhsSource` now live only under `clode.problem._core`; the curated public problem API centers `InitialValueProblem` and the authoring/conversion helpers instead.
- Built-in observer definitions, resolved observer specs, and observer-state naming are now aligned across the Python host layer, OpenCL metadata, and the active kernels.
- Execution-setting defaults and compatibility resolution now flow through one canonical solver-settings path, and simulators keep internal copies of caller-provided `SolverParams` bundles instead of aliasing them.
- Built-in stepper definitions, traits, and OpenCL build mapping now resolve through a Python-owned stepper-definition catalog instead of raw registry tables.
- The current package still has a few cleanup targets, especially shared kernel-math helpers and component tests for observer/stepper internals, continuation-state semantics for diverged work-item times, and later chunking work on top of the landed IVP, solver-state, output-policy, observer, execution-setting, and stepper-definition boundaries.

## Session-Start Guidance

- Prefer canonical imports from `clode`, `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- Treat root modules such as `clode.solver`, `clode.features`, `clode.trajectory`, `clode.types`, `clode.function_converter`, `clode.xpp_parser`, and `clode.opencl_builtins` as compatibility barrels only.
- New implementation work should land in the canonical packages, not in the compatibility barrels.
- Use `.design/next_pr.md` to decide what to work on next; use this file to understand where code now belongs.

## Top-Level Repo Map

| Path | Role | Notes |
| --- | --- | --- |
| `clode/` | Maintained package source | All runtime-critical code lives here |
| `docs/` | MkDocs source | Current user-facing Python docs |
| `examples/` | Current example scripts and sample RHS files | Live Python/OpenCL examples |
| `test/` | Supported regression suite | Organized by bundle and marker rather than final directory layout |
| `tools/` | Development and diagnostics helpers | Includes the test-bundle runner and OpenCL probe |
| `.design/archived/` | Archived migration, cleanup, and bug notes | High-signal design history, but some statements are historical only |
| `dist/`, `site/`, `clode.egg-info/` | Generated packaging or docs artifacts | Useful locally, not authoritative source |
| `paper/` | Paper describing this project (target: Journal of Open Source Software) | Orthogonal to the runtime and package internals |

## Package Map

### Public layer

- `clode/__init__.py`: curated export surface and package version.
- `clode/problem/*`: user-facing problem definition via `InitialValueProblem` plus Python/OpenCL/XPP authoring helpers. Lower-level support types such as `ProblemInfo` and `RhsSource` live only in the internal `_core` module as derived support concepts.
- `clode/runtime/*`: public OpenCL query, device selection, and stdlib-logging helpers.
- `clode/simulation/*`: solver params, simulator orchestration, current ensemble helpers, and result containers.
- `clode/observers/*`: observer enums, parameter schema, public feature-name helpers, and the internal built-in observer-definition catalog in `_definitions.py`.
- Root compatibility barrels remain for public flat imports such as `solver.py`, `features.py`, `trajectory.py`, `function_converter.py`, `xpp_parser.py`, and `opencl_builtins.py`.

### Internal OpenCL implementation

- `clode/_opencl/models.py`: build keys, source bundles, problem shape, and precision enums.
- `clode/_opencl/runtime.py`: explicit single-device context and queue creation plus `RuntimeSelection` normalization into a concrete PyOpenCL runtime.
- `clode/_opencl/registry.py`: entrypoint mapping plus stepper and observer resolution routed through semantic definition catalogs.
- `clode/_opencl/source_builder.py`: kernel assembly, build options, and kernel-tree digesting.
- `clode/_opencl/program_cache.py`: runtime-scoped OpenCL program cache.
- `clode/_opencl/structs.py`: device-matched struct dtypes.
- `clode/_opencl/buffers.py`: buffer allocation plus flatten/reshape rules.
- `clode/_opencl/observer_metadata.py`: runtime-specific matched observer-state struct resolution driven by the observer-definition catalog.
- `clode/_opencl/executors.py`: transient, trajectory, and feature executors constructed directly by the simulator layer.

The historical backend-shim, protocol/factory facade, and binding-named compatibility layers have been removed.

## Compatibility Surface

- `clode.__init__` remains the main stable top-level import surface and is worth preserving.
- The flat root files are packaging-neutral compatibility shims. They are not required by setuptools package discovery, package data shipping, or wheel/sdist correctness.
- Their real value is import-path continuity for downstream users, examples, old notes, type annotations, any serialized or pickled objects that still mention historical module paths, and collaborator orientation while the semantic package layout continues to settle.
- Their main cost is duplicate API surface, extra documentation burden, and a greater chance that future work accidentally lands in the wrong module.
- If they are removed later, do it as a normal deprecation cycle: shift docs/examples/tests first, optionally add warnings, then remove them in a deliberate release.

### Kernel tree

- `clode/kernels/transient.cl`: base transient kernel.
- `clode/kernels/trajectory.cl`: trajectory storage path.
- `clode/kernels/initializeObserver.cl` and `clode/kernels/features.cl`: feature and observer lifecycle.
- `clode/kernels/steppers/*.clh`: explicit and adaptive stepper implementations.
- `clode/kernels/observers/*.clh`: built-in observer implementations.
- `clode/kernels/odedriver.cl`: currently unused unified-driver prototype; intentionally kept as deferred future-design context, not active runtime code.

## Runtime Model

- One OpenCL work-item advances one ODE instance.
- Problem arrays are flattened in Fortran order on the Python side.
- Kernel builds remain compile-time specialized by precision, stepper, observer, and problem dimensions.
- Simulators now orchestrate an `InitialValueProblem`, runtime binding, continuation state, output policy, and compatibility helper methods such as simulator-side `set_ensemble()` delegation.
- `FeatureSimulator` uses persistent observer state plus an optional warmup kernel for two-pass observers.
- The solver now owns the authoritative continuation time base. Observer finalizers no longer rebase time, and exact split-window continuation is defined by continuing from the attained per-item `tf`.
- A first-pass internal solver-state boundary now lives in `clode/simulation/_state.py`, while `_opencl/executors.py` treats host-side mirrors as transfer caches rather than semantic owners. Current runtime continuation still uses a shared requested `tspan` plus per-item `dt`, attained `tf`, and RNG continuation rather than a full device-side per-work-item state object.
- Exact absolute-time continuation still requires caller-managed `t_span`. For fixed-step runs, the robust handoff point is the attained `tf` from `get_final_time()`, not the requested endpoint.

## Current Constraints And Live Debt

- The runtime is explicitly single-device only. Any future multi-device execution would require a dedicated API and execution model rather than reviving removed transition-era selectors.
- The first semantic IVP pass is now live, and the curated public problem API now centers `InitialValueProblem`; lower-level helper types such as `ProblemInfo` and `RhsSource` remain internal support concepts under `clode.problem._core`.
- Simulators still carry compatibility delegates for batch reshaping and maintain cached problem arrays alongside the IVP; richer batch-generation helpers (`grid`, random, quasi-random) have not yet been added around the IVP model.
- Solver-related execution state now has a clearer internal home in `clode/simulation/_state.py` and executor transfer caches, and integration policy is now explicitly split from trajectory output/storage policy through internal settings views, OpenCL buffers, and kernel-facing structs. Public `SolverParams` remains a compatibility bundle at the simulator boundary.
- Common continuation state now has a Python-owned first pass via `SolverState`, but there is still no device-side per-work-item `t0` or richer completion/error status model.
- Exact shared-window continuation is only representable when the ensemble agrees on one attained `tf`; with diverged per-work-item final times, any next shared `t_span` is an explicit approximation or policy choice rather than one exact continuation update.
- `SolverParams` still mixes integration controls with trajectory-storage controls (`max_store`, `nout`) in the public compatibility bundle even though the internal and kernel-facing execution path now treats them separately.
- `FeatureSimulator` still exposes a broad legacy `observer_*` scalar compatibility surface alongside `ObserverParams`, even though the observer-definition model is now the clearer internal semantic boundary.
- `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` now resolve scalar solver arguments and prebuilt `SolverParams` bundles through the same canonical helper, but `SolverParams` still mixes integration controls with trajectory-output policy at the public compatibility boundary.
- Built-in stepper definitions now live in `clode/simulation/_stepper_definitions.py`, and runtime validation plus source-building now consume those definitions instead of raw stepper string tables.
- `ObserverParams` now owns the canonical built-in observer defaults, but the legacy constructor alias names in `FeatureSimulator` remain as the compatibility surface.
- Built-in observer definitions now live in `clode/observers/_definitions.py`, and shared resolved observer specs now drive feature-build defines, build-key selection, and feature-executor invalidation boundaries.
- Optional event storage still participates in compile-time observer layout and buffer sizing, but runtime observer settings now exclude event timestamp capacity and treat it as explicit output/layout policy instead.
- `shift_tspan()` still advances the requested window rather than the attained final time, so exact absolute-time continuation remains a caller-managed policy.
- Shared OpenCL numerical helpers for compensated accumulation and structured time updates are still mostly TODOs or local kernel code rather than one small reusable utility layer.
- Fixed-step kernels still advance time with `ti += dt`, so very long absolute-time runs with small `dt` remain precision-sensitive; evaluating `t0 + step * dt` or related structured-time models is still future work.
- RNG continuation details are persisted in separate common buffers rather than a clearer per-work-item state object, which will matter again when evaluating Random123.
- `.design/archived/` is useful for rationale and bug archaeology, but some archived statements about the old wrapper path are now historical only.

## Test And Tooling Quick Lookup

- `tools/run_test_bundle.py`: authoritative bundle map.
- `test/core_numerics/`: exact-solution and kernel-level regression backbone.
- `test/core_numerics/test_stochastic.py`: seeded stochastic continuation and Ornstein-Uhlenbeck stationary-moment coverage.
- `test/core_numerics/test_features_basicall.py`: `basicall` exact-statistics and split-window continuation coverage.
- `test/test_simulation_contracts.py`: current simulation and observer behavior contracts.
- There is not yet a dedicated kernel-component test layer between the end-to-end numerical regressions and the `_opencl` support-layer tests.
- `test/test_problem_rhs_source.py`: RHS source-ingestion and digest coverage.
- `test/test_opencl_models.py`, `test/test_opencl_source_builder.py`, `test/test_opencl_runtime.py`, `test/test_opencl_buffers.py`, `test/test_opencl_structs.py`: canonical internal OpenCL support-layer tests.
- `tools/probe_opencl_runtime.py`: distinguishes runtime/compiler failures from clODE kernel failures.
- `docs/init_runtime.md`: current single-device runtime-selection story.

## Packaging Quick Lookup

- `pyproject.toml`: authoritative packaging and dependency configuration.
- `MANIFEST.in`: current sdist include and prune rules.
- Kernel assets ship as package data.
- Docs and tests are not required for runtime execution. Tests are still useful for downstream verification, while docs are the easier thing to omit if sdist slimming becomes desirable.
- From a packaging perspective, the flat compatibility barrels are optional; they exist only to preserve import compatibility, not because the build needs them.
