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

- [ ] P2 Device-side per-work-item current time / `t0` if exact continuation after diverged `tf` becomes a real priority. depends: continuation-state semantics and guardrails. refs: `clode/simulation/_state.py`, `clode/_opencl/executors.py`, `clode/_opencl/buffers.py`, `clode/kernels/transient.cl`, `.design/reference/continuation_timebase_note.md`
- [ ] P1 Batch-generation helpers such as `grid`, random, and quasi-random sampling layered on top of the IVP model and current shape metadata instead of keeping ensemble creation buried in `Simulator`. refs: `clode/simulation/base.py`, `.design/reference/semantic_layout_audit.md`
- [ ] P2 Execution-model experiments: `per-work-item` vs `per-work-group` vs shared/global solver state. refs: `clode/kernels/odedriver.cl`

## Observer And Feature Model

- [ ] P2 Add a small set of dynamical-systems-oriented observers or features such as direction-of-crossing or Poincare-section style events once the built-in observer model settles further. refs: `clode/observers/_definitions.py`, `clode/kernels/observers/`, `docs/examples.md`
- [ ] P2 Observer-specific parameter models/classes instead of one broad `ObserverParams`. depends: explicit observer-definition model. refs: `clode/features.py`, `clode/kernels/observers.cl`
- [ ] P2 Support aux variables as event/feature variables. depends: explicit observer-definition model. refs: `clode/kernels/observers.cl`
- [ ] P2 Custom/composable observers from Python-authored definitions. depends: explicit observer-definition model and codegen story. refs: `clode/features.py`, `clode/function_converter.py`

## Numerical Methods And Kernel Math

- [ ] P0 Empirical single-precision numerics demonstrations and docs before broader mitigation rollout: compare `runningMeanTime(...)` against compensated integral accumulation, compare direct time addition against compensated or structured time updates, compare sampled threshold timestamps against interpolation alternatives, and capture the results in public docs and helper tests. refs: `examples/single_precision_accuracy.py`, `docs/numerical_accuracy.md`, `.design/reference/single_precision_numerics_note.md`, `clode/kernels/clODE_utilities.cl`, `test/kernel_components/test_kernel_math.py`, `clode/kernels/steppers/`
- [ ] P1 Broaden shared numerical-helper adoption beyond the initial helper foundation: carry compensated mean bookkeeping into the remaining observers and resume stepper time-base work only where the empirical demos justify the extra complexity, preferring a dual-realtype compensated time path plus relative elapsed bookkeeping over more single-float variants. depends: empirical single-precision numerics demonstrations and docs. refs: `clode/kernels/clODE_utilities.cl`, `clode/kernels/observers/`, `clode/kernels/steppers/`, `test/kernel_components/test_kernel_math.py`
- [ ] P1 Better trajectory/output modes: chunking, variable subsets, specified output times, dense output. depends: integration/output separation. refs: `clode/trajectory.py`, `clode/kernels/trajectory.cl`, `.design/reference/chunked_execution_audit.md`
- [ ] P2 Time-base and interpolation accuracy helpers (dual-realtype time bookkeeping, `TwoSum`, `t0 + step * dt`, step-counter width and overflow budget, inverse-linear and slope-aware threshold timestamps, three-sample local-extremum helpers, dense interpolants, fixed-step endpoint handling). depends: explicit per-work-item `t0`. refs: `clode/kernels/clODE_utilities.cl`, `clode/kernels/realtype.cl`, `.design/reference/single_precision_numerics_note.md`, `.design/archived/pre_backend_readiness_2026_05_05/fixed_step_endpoint_bug_audit.md`
- [ ] P3 Rounding-mode audit for low-precision kernels: decide whether round-to-nearest assumptions are sufficient and whether stochastic-rounding experiments are worthwhile given reproducibility and OpenCL support constraints. depends: empirical single-precision numerics demonstrations and docs. refs: `.design/reference/single_precision_numerics_note.md`, `clode/kernels/clODE_utilities.cl`, `clode/kernels/steppers/`
- [ ] P2 Better adaptive-step controllers (I/PI/PID). refs: `clode/kernels/steppers.cl`
- [ ] P2 Bounds/blow-up detection in solver structs/kernels. refs: `clode/kernels/clODE_struct_defs.cl`
- [ ] P3 Implicit fixed/adaptive steppers for stiff systems. depends: explicit solver-state model and likely Jacobian story. refs: `clode/kernels/steppers.cl`, `.design/development_roadmap.md`
- [ ] P3 Jacobian generation/emission for RHS. depends: stronger RHS IR. blocks: robust implicit solver work. refs: `.design/ideas.md`, `clode/function_converter.py`

## OpenCL Program And Runtime Model

- [ ] P1 Solution-buffer / solver-state kernel abstraction shared by steppers and observers. depends: explicit solver-state model. refs: `clode/kernels/observers.cl`, `clode/kernels/odedriver.cl`
- [ ] P1 Kernel specialization and source-assembly audit: decide whether `KernelKind`, entrypoint selection, and the current `#define`/`#include` model should specialize more aggressively, but keep kernel files separate from Python definition catalogs unless stronger evidence appears. refs: `clode/_opencl/source_builder.py`, `clode/_opencl/registry.py`, `clode/_opencl/models.py`, `clode/kernels/`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P1 Leverage remaining PyOpenCL runtime/build helpers (`cache_dir`, broader `characterize` helpers, `capture_call`) plus device-side fill/map helpers before growing more custom diagnostics or transfer code. refs: `clode/runtime/query.py`, `clode/_opencl/runtime.py`, `clode/_opencl/program_cache.py`, `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Benchmark whether PyOpenCL `MemoryPool` / `ImmediateAllocator` reduce transient, trajectory, or feature buffer churn before inventing custom allocation policy. depends: clearer PyOpenCL leverage strategy. refs: `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Audit whether PyOpenCL array-side helpers (`Array`, `ElementwiseKernel`, reductions, scans) belong only in diagnostics, preprocessing, or testing helpers instead of the solver hot path. depends: clearer PyOpenCL leverage strategy. refs: `clode/_opencl/buffers.py`, `clode/_opencl/executors.py`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Audit underused OpenCL kernel features (for example vector types) for measurable wins in clODE kernels without obscuring solver semantics. refs: `clode/kernels/`, `.design/development_roadmap.md`
- [ ] P2 Compatibility-barrel deprecation plan: decide whether to keep only `clode.__init__` as the stable top-level barrel and retire flat modules such as `clode.solver`, `clode.features`, `clode.trajectory`, `clode.types`, `clode.function_converter`, `clode.xpp_parser`, and `clode.opencl_builtins`. depends: docs/examples shifted to canonical imports and deprecation appetite agreed. refs: `.design/package_state.md`, `.design/development_roadmap.md`, `pyproject.toml`
- [ ] P3 Dedicated multi-device execution design only if concrete demand appears; do not overload the current single-device runtime selectors. refs: `docs/init_runtime.md`, `clode/runtime/selection.py`, `.design/development_roadmap.md`

## RHS IR, Conversion, And Interop

- [ ] P2 Stronger internal RHS representation and simplification passes. blocks: Jacobian emission, round-tripping, better interop. refs: `clode/function_converter.py`, `.design/ideas.md`
- [ ] P2 Generic RHS interop with SciPy/other solver packages and fewer converter restrictions across Python, OpenCL, and XPP-defined problems. depends: stronger RHS IR. refs: `clode/problem/ivp.py`, `clode/function_converter.py`, `docs/specifying_odes.md`
- [ ] P2 Finish the lark-based XPP parser path and configuration handling. refs: `clode/xpp_parser.py`
- [ ] P3 Export paths from internal RHS form back to Python/XPP. depends: stronger RHS IR. refs: `clode/xpp_parser.py`, `clode/function_converter.py`

## Performance, Randomness, And Scaling

- [ ] P1 Ensemble batching for device-capacity limits. depends: IVP-owned batch semantics and integration/output separation. refs: `clode/trajectory.py`, `.design/development_roadmap.md`, `.design/reference/chunked_execution_audit.md`
- [ ] P2 External RNG library audit (for example PyOpenCL `clrandom` or kernel-side Random123 headers) plus whether the current Box-Muller/prepared-Wiener continuation fields should become a named per-work-item RNG state. depends: explicit solver-state model and reproducibility requirements. refs: `clode/_opencl/executors.py`, `clode/_opencl/buffers.py`, `clode/kernels/clODE_random.cl`, `clode/kernels/transient.cl`, `clode/kernels/features.cl`, `.design/reference/pyopencl_leverage_audit.md`
- [ ] P2 Write-avoidance and memory-traffic audit for trajectory/event storage. depends: chunked/output redesign. refs: `clode/kernels/trajectory.cl`, `clode/kernels/odedriver.cl`

## Testing, Diagnostics, And Scope

- [ ] P1 Expand kernel-component coverage beyond helper math and the current basic/basicall observer contracts into build-key invalidation, observer storage, and stepper execution contracts. refs: `test/kernel_components/`, `clode/_opencl/source_builder.py`, `.design/reference/testing_audit.md`
- [ ] P1 Cache/build-key/invalidation coverage around program rebuilds and observer changes. refs: `clode/_opencl/program_cache.py`, `clode/_opencl/source_builder.py`, `clode/features.py`
- [ ] P1 Test-surface audit for drift, weak success conditions, and unnecessary coverage before more large internal refactors. refs: `test/`, `tools/run_test_bundle.py`, `.design/reference/testing_audit.md`
- [ ] P2 Package/API hygiene: typing cleanup, stepper-specific parameter subsets, and trajectory output ergonomics once continuation policy and stepper-specific public semantics are clearer. refs: `clode/solver.py`, `clode/features.py`, `clode/trajectory.py`
- [ ] P2 Examples or docs that show where compensated summation and related single-precision safeguards help in practice. refs: `examples/`, `docs/examples.md`, `docs/performance_notes.md`, `clode/kernels/clODE_utilities.cl`
- [ ] P1 Runtime-install guidance: reuse upstream PyOpenCL install/runtime guidance where it reduces local duplication and point users to optional runtime packaging where appropriate. refs: `docs/install.md`, `docs/querying_opencl.md`, `.design/reference/public_surfaces_plan.md`, `.design/reference/docs_layout_plan.md`
- [ ] P2 Public-facing package narrative and publication-readiness: maintain the shared landing page, keep the docs IA iterative (`Getting Started` / `Guides` / narrow `Examples` / `Reference`) while API changes settle, and defer citation metadata plus release-tag hygiene until the scientific, numerical, performance, runtime, and UX surfaces are stronger. refs: `README.md`, `CONTRIBUTING.md`, `docs/index.md`, `docs/examples.md`, `docs/performance_notes.md`, `paper/paper.md`, `.design/reference/project_principles.md`, `.design/reference/joss_audit.md`, `.design/reference/public_surfaces_plan.md`, `.design/reference/docs_layout_plan.md`
- [ ] P2 Scope audit: what belongs in clODE vs sibling/helper packages. refs: `.design/package_state.md`, `.design/development_roadmap.md`, `.design/reference/project_principles.md`

## Inbox

Temporary holding area for rough notes. Process this section with the `design-ideas-inbox` skill when routing items into the maintained `.design` framework.
