# clODE Planning Board

Purpose: terse backlog and sequencing board.
Read when: you need dependencies, adjacent work, or follow-on options beyond the active PR target.
Update when: work is added, reprioritized, completed, split, or deduplicated.

Keep this file terse and editable.

Companions:

- `.design/package_state.md`: factual repo and package map
- `.design/development_roadmap.md`: longer rationale and prioritization
- `.design/next_pr.md`: current narrow implementation target
- `.design/reference/README.md`: focused deep dives and audit notes
- `.design/archived/`: historical detail and design archaeology

Format:

- `[ ] P0/P1/P2 item. depends: prerequisite. blocks: downstream work. refs: code/doc paths`
- keep details out of this file when a single line is enough
- prefer updating an existing line over adding a near-duplicate

## Core Execution And State Semantics

- [x] P0 InitialValueProblem-first batch semantics at the simulator boundary: defaults, batch shaping, remembered result shape, simulator `ivp=` construction, and Python-backed SciPy-style callability now live around `InitialValueProblem` so simulators read more clearly as orchestration objects. refs: `.design/next_pr.md`, `clode/problem/ivp.py`, `clode/simulation/base.py`, `test/test_initial_value_problem.py`, `test/test_initial_value_problem_runtime.py`
- [x] P1 IVP public-surface cleanup: `InitialValueProblem` is now the promoted user-facing problem API in the curated docs and exports, while `ProblemInfo`/`RhsSource` live only in the internal `clode.problem._core` support layer. refs: `docs/api_reference.md`, `docs/specifying_odes.md`, `clode/problem/__init__.py`, `clode/problem/_core.py`, `.design/reference/ivp_api_test_plan.md`
- [x] P0 First-pass explicit solver-state and cache ownership cleanup: IVP now owns next-solve problem data, `clode/simulation/_state.py` owns Python-side solver state and fetched-output caches, and `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners. refs: `.design/reference/continuation_timebase_note.md`, `.design/reference/solver_state_implementation_plan.md`, `clode/simulation/_state.py`, `clode/_opencl/executors.py`, `test/test_opencl_executors.py`, `test/test_simulation_contracts.py`
- [ ] P1 Public continuation-policy helper/API (`requested-window` vs attained-`tf` continuation) so callers do not have to hand-roll `set_tspan(get_final_time())`. depends: observer-definition cleanup and stepper-definition cleanup. refs: `docs/continuation.md`, `clode/simulation/base.py`
- [x] P0 Separate integration state from output/storage policy (`max_store`, `nout`, event storage) so trajectory and feature allocation stop reading like solver state. refs: `clode/simulation/params.py`, `clode/observers/types.py`, `clode/_opencl/buffers.py`, `clode/_opencl/structs.py`, `clode/_opencl/executors.py`, `clode/kernels/clODE_struct_defs.cl`, `test/test_opencl_executors.py`, `test/test_opencl_buffers.py`, `test/test_opencl_structs.py`
- [ ] P1 Python-owned stepper-definition model: explicit fixed/adaptive, explicit/implicit, deterministic/stochastic traits plus OpenCL mapping instead of splitting stepper semantics across `Stepper`, `_opencl.registry`, and kernel defines. depends: observer-definition cleanup. blocks: cleaner implicit-stepper work, more intuitive source assembly, and later public config redesign. refs: `clode/simulation/base.py`, `clode/_opencl/registry.py`, `clode/kernels/steppers.cl`, `.design/reference/semantic_layout_audit.md`
- [ ] P1 Batch-generation helpers such as `grid`, random, and quasi-random sampling layered on top of the IVP model and current shape metadata instead of keeping ensemble creation buried in `Simulator`. refs: `clode/simulation/base.py`, `.design/reference/semantic_layout_audit.md`
- [ ] P2 Execution-model experiments: `per-work-item` vs `per-work-group` vs shared/global solver state. refs: `clode/kernels/odedriver.cl`

## Observer And Feature Model

- [ ] P0 Explicit observer-definition model: params, persistent state, event storage, warmup, feature names. blocks: custom/composable observers, observer-state cleanup, and later public config redesign. refs: `.design/next_pr.md`, `clode/observers/metadata.py`, `clode/_opencl/observer_metadata.py`, `clode/_opencl/executors.py`, `clode/kernels/observers.cl`
- [ ] P1 Separate persistent observer state from optional event-output capacity. depends: explicit observer-definition model. refs: `clode/_opencl/observer_metadata.py`, `clode/kernels/features.cl`, `.design/development_roadmap.md`
- [ ] P2 Observer-specific parameter models/classes instead of one broad `ObserverParams`. depends: explicit observer-definition model. refs: `clode/features.py`, `clode/kernels/observers.cl`
- [ ] P2 Support aux variables as event/feature variables. depends: explicit observer-definition model. refs: `clode/kernels/observers.cl`
- [ ] P2 Custom/composable observers from Python-authored definitions. depends: explicit observer-definition model and codegen story. refs: `clode/features.py`, `clode/function_converter.py`

## Numerical Methods And Kernel Math

- [ ] P1 Better trajectory/output modes: chunking, variable subsets, specified output times, dense output. depends: integration/output separation. refs: `clode/trajectory.py`, `clode/kernels/trajectory.cl`, `.design/reference/chunked_execution_audit.md`
- [ ] P2 Time-base and interpolation accuracy helpers (`TwoSum`, `t0 + step * dt`, dense interpolants, fixed-step endpoint handling). depends: explicit per-work-item `t0`. refs: `clode/kernels/clODE_utilities.cl`, `clode/kernels/realtype.cl`, `.design/archived/pre_backend_readiness_2026_05_05/fixed_step_endpoint_bug_audit.md`
- [ ] P2 Better adaptive-step controllers (I/PI/PID). refs: `clode/kernels/steppers.cl`
- [ ] P2 Bounds/blow-up detection in solver structs/kernels. refs: `clode/kernels/clODE_struct_defs.cl`
- [ ] P3 Implicit fixed/adaptive steppers for stiff systems. depends: explicit solver-state model and likely Jacobian story. refs: `clode/kernels/steppers.cl`, `.design/development_roadmap.md`
- [ ] P3 Jacobian generation/emission for RHS. depends: stronger RHS IR. blocks: robust implicit solver work. refs: `.design/ideas.md`, `clode/function_converter.py`

## OpenCL Program And Runtime Model

- [ ] P1 Solution-buffer / solver-state kernel abstraction shared by steppers and observers. depends: explicit solver-state model. refs: `clode/kernels/observers.cl`, `clode/kernels/odedriver.cl`
- [ ] P2 Audit Python-level source assembly vs current `#define`/`#include` model. refs: `clode/_opencl/source_builder.py`, `clode/kernels/`, `.design/archived/pyopencl_cleanup_closeout_2026_05_08/pyopencl_backend_design.md`
- [ ] P1 Leverage remaining PyOpenCL runtime/build helpers (`cache_dir`, broader `characterize` helpers, `capture_call`) plus device-side fill/map helpers before growing more custom diagnostics or transfer code. refs: `clode/runtime/query.py`, `clode/_opencl/runtime.py`, `clode/_opencl/program_cache.py`, `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Benchmark whether PyOpenCL `MemoryPool` / `ImmediateAllocator` reduce transient, trajectory, or feature buffer churn before inventing custom allocation policy. depends: clearer PyOpenCL leverage strategy. refs: `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Audit whether PyOpenCL array-side helpers (`Array`, `ElementwiseKernel`, reductions, scans) belong only in diagnostics, preprocessing, or testing helpers instead of the solver hot path. depends: clearer PyOpenCL leverage strategy. refs: `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Audit underused OpenCL kernel features (for example vector types) for measurable wins in clODE kernels without obscuring solver semantics. refs: `clode/kernels/`, `.design/development_roadmap.md`
- [ ] P2 Compatibility-barrel deprecation plan: decide whether to keep only `clode.__init__` as the stable top-level barrel and retire flat modules such as `clode.solver`, `clode.features`, `clode.trajectory`, `clode.types`, `clode.function_converter`, `clode.xpp_parser`, and `clode.opencl_builtins`. depends: docs/examples shifted to canonical imports and deprecation appetite agreed. refs: `.design/package_state.md`, `.design/development_roadmap.md`, `pyproject.toml`
- [ ] P3 Dedicated multi-device execution design only if concrete demand appears; do not overload the current single-device runtime selectors. refs: `docs/init_runtime.md`, `clode/runtime/selection.py`, `.design/development_roadmap.md`

## RHS IR, Conversion, And Interop

- [ ] P2 Stronger internal RHS representation and simplification passes. blocks: Jacobian emission, round-tripping, better interop. refs: `clode/function_converter.py`, `.design/ideas.md`
- [ ] P2 Generic RHS interop with SciPy/other solver packages and fewer converter restrictions across Python, OpenCL, and XPP-defined problems. depends: stronger RHS IR. refs: `clode/function_converter.py`, `.design/next_pr.md`
- [ ] P2 Finish the lark-based XPP parser path and configuration handling. refs: `clode/xpp_parser.py`
- [ ] P3 Export paths from internal RHS form back to Python/XPP. depends: stronger RHS IR. refs: `clode/xpp_parser.py`, `clode/function_converter.py`

## Performance, Randomness, And Scaling

- [ ] P1 Ensemble batching for device-capacity limits. depends: IVP-owned batch semantics and integration/output separation. refs: `clode/trajectory.py`, `.design/development_roadmap.md`, `.design/reference/chunked_execution_audit.md`
- [ ] P2 External RNG library audit (for example PyOpenCL `clrandom` or kernel-side Random123 headers) plus whether the current Box-Muller/prepared-Wiener continuation fields should become a named per-work-item RNG state. depends: explicit solver-state model and reproducibility requirements. refs: `clode/_opencl/executors.py`, `clode/_opencl/buffers.py`, `clode/kernels/clODE_random.cl`, `clode/kernels/transient.cl`, `clode/kernels/features.cl`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Write-avoidance and memory-traffic audit for trajectory/event storage. depends: chunked/output redesign. refs: `clode/kernels/trajectory.cl`, `clode/kernels/odedriver.cl`

## Testing, Diagnostics, And Scope

- [x] P1 Logging/runtime diagnostics audit: remove the log-level compatibility API, use standard logging with a minimal `configure_logging(...)` / `get_logger(...)` surface, keep explicit `print_*` helpers explicit, and lean on PyOpenCL diagnostics instead of growing custom logging glue. refs: `.design/archived/post_logging_cleanup_2026_05_12/logging_audit.md`, `clode/runtime/logging.py`, `clode/runtime/query.py`, `clode/_opencl/runtime.py`, `clode/_opencl/program_cache.py`
- [ ] P1 Kernel-component tests with tiny synthetic models and kernels beyond end-to-end simulator tests. refs: `test/core_numerics/`, `clode/_opencl/source_builder.py`, `.design/reference/testing_audit.md`
- [ ] P1 Cache/build-key/invalidation coverage around program rebuilds and observer changes. refs: `clode/_opencl/program_cache.py`, `clode/_opencl/source_builder.py`, `clode/features.py`
- [ ] P2 Package/API hygiene: centralize defaults, typing cleanup, stepper-specific parameter subsets, trajectory output ergonomics. refs: `clode/solver.py`, `clode/features.py`, `clode/trajectory.py`
- [ ] P1 Runtime-install guidance: reuse upstream PyOpenCL install/runtime guidance where it reduces local duplication and point users to optional runtime packaging where appropriate. refs: `docs/install.md`, `docs/querying_opencl.md`, `.design/reference/public_surfaces_plan.md`, `.design/reference/docs_layout_plan.md`
- [ ] P1 Public-facing package narrative and publication-readiness: maintain the shared landing page, keep the docs IA iterative (`Getting Started` / `Guides` / narrow `Examples` / `Reference`) while API changes settle, add citation and contributor metadata, and keep the paper aligned with current supported workflows. refs: `README.md`, `CONTRIBUTING.md`, `docs/index.md`, `docs/examples.md`, `docs/performance_notes.md`, `paper/paper.md`, `.design/reference/project_principles.md`, `.design/reference/joss_audit.md`, `.design/reference/public_surfaces_plan.md`, `.design/reference/docs_layout_plan.md`
- [ ] P2 Scope audit: what belongs in clODE vs sibling/helper packages. refs: `.design/package_state.md`, `.design/development_roadmap.md`, `.design/reference/project_principles.md`

## Inbox

Temporary holding area for rough notes. Process this section with the `design-ideas-inbox` skill when routing items into the maintained `.design` framework.

- none currently
