# clODE Planning Board

Short living board. Keep this file terse and editable.

Companions:

- `tmp/package_state.md`: factual repo and package map
- `tmp/development_roadmap.md`: longer rationale and prioritization
- `tmp/next_pr.md`: current narrow implementation target
- `tmp/archived/`: historical detail and design archaeology

Format:

- `[ ] P0/P1/P2 item. depends: prerequisite. blocks: downstream work. refs: code/doc paths`
- keep details out of this file when a single line is enough
- prefer updating an existing line over adding a near-duplicate

## Core Execution And State Semantics

- [ ] P0 Continuation semantics and explicit solver/observer state. blocks: chunking, batching, implicit methods, non-autonomous correctness. refs: `tmp/next_pr.md`, `clode/solver.py`, `clode/kernels/features.cl`, `tmp/archived/pyopencl_cleanup_closeout_2026_05_08/python_simulation_flow_reference.md`
- [ ] P1 Separate integration state from output/storage policy (`max_store`, `nout`, event storage). depends: continuation semantics. blocks: chunked trajectory storage, observer-storage cleanup. refs: `clode/trajectory.py`, `clode/kernels/trajectory.cl`, `tmp/development_roadmap.md`
- [ ] P1 Ensemble resizing/broadcast semantics for `set_ensemble()` and `set_repeat_ensemble()`. depends: clearer solver state model. refs: `clode/solver.py`
- [ ] P2 Execution-model experiments: `per-work-item` vs `per-work-group` vs shared/global solver state. depends: clearer solver-state model. refs: `clode/kernels/odedriver.cl`

## Observer And Feature Model

- [ ] P1 Explicit observer-definition model: params, persistent state, event storage, warmup, feature names. depends: continuation semantics. blocks: custom/composable observers. refs: `clode/features.py`, `clode/observers/metadata.py`, `clode/_opencl/observer_metadata.py`, `clode/kernels/observers.cl`
- [ ] P1 Separate persistent observer state from optional event-output capacity. depends: explicit observer-definition model. refs: `clode/_opencl/observer_metadata.py`, `clode/kernels/features.cl`, `tmp/development_roadmap.md`
- [ ] P2 Observer-specific parameter models/classes instead of one broad `ObserverParams`. depends: explicit observer-definition model. refs: `clode/features.py`, `clode/kernels/observers.cl`
- [ ] P2 Support aux variables as event/feature variables. depends: explicit observer-definition model. refs: `clode/kernels/observers.cl`
- [ ] P2 Custom/composable observers from Python-authored definitions. depends: explicit observer-definition model and codegen story. refs: `clode/features.py`, `clode/function_converter.py`

## Numerical Methods And Kernel Math

- [ ] P1 Better trajectory/output modes: chunking, variable subsets, specified output times, dense output. depends: integration/output separation. refs: `clode/trajectory.py`, `clode/kernels/trajectory.cl`
- [ ] P2 Time-accumulation and interpolation accuracy helpers (`TwoSum`, dense interpolants, fixed-step endpoint handling). refs: `clode/kernels/clODE_utilities.cl`, `clode/kernels/realtype.cl`, `tmp/archived/pre_backend_readiness_2026_05_05/fixed_step_endpoint_bug_audit.md`
- [ ] P2 Better adaptive-step controllers (I/PI/PID). refs: `clode/kernels/steppers.cl`
- [ ] P2 Bounds/blow-up detection in solver structs/kernels. refs: `clode/kernels/clODE_struct_defs.cl`
- [ ] P3 Implicit fixed/adaptive steppers for stiff systems. depends: continuation semantics, solver-state model, likely Jacobian story. refs: `clode/kernels/steppers.cl`, `tmp/development_roadmap.md`
- [ ] P3 Jacobian generation/emission for RHS. depends: stronger RHS IR. blocks: robust implicit solver work. refs: `tmp/ideas.md`, `clode/function_converter.py`

## OpenCL Program And Runtime Model

- [ ] P1 Solution-buffer / solver-state kernel abstraction shared by steppers and observers. depends: continuation semantics. refs: `clode/kernels/observers.cl`, `clode/kernels/odedriver.cl`
- [ ] P2 Audit Python-level source assembly vs current `#define`/`#include` model. refs: `clode/_opencl/source_builder.py`, `clode/kernels/`, `tmp/archived/pyopencl_cleanup_closeout_2026_05_08/pyopencl_backend_design.md`
- [ ] P1 Leverage PyOpenCL compiler cache, `MemoryPool`, `capture_call`, and `characterize` utilities before adding more custom runtime helpers. refs: `clode/_opencl/runtime.py`, `clode/_opencl/program_cache.py`, `clode/_opencl/buffers.py`, `tmp/pyopencl_leverage_audit.md`
- [ ] P1 Role-based package layout migration after the PyOpenCL-first decision: collapse the transition-era `_backends/` seam, keep root compatibility re-exports, and reorganize around `problem/`, `runtime/`, `simulation/`, `observers/`, and `_opencl/`. depends: consensus on `tmp/module_layout_plan.md`. refs: `tmp/module_layout_plan.md`, `tmp/backend_strategy_audit.md`, `tmp/package_state.md`
- [ ] P3 Decide whether `device_ids` becomes real multi-device work or is retired. depends: explicit multi-device scope decision. refs: `docs/init_runtime.md`, `clode/runtime.py`

## RHS IR, Conversion, And Interop

- [ ] P2 Stronger internal RHS representation and simplification passes. blocks: Jacobian emission, round-tripping, better interop. refs: `clode/function_converter.py`, `tmp/ideas.md`
- [ ] P2 Python RHS interop with SciPy/other solver packages and fewer converter restrictions. depends: stronger RHS IR. refs: `clode/function_converter.py`
- [ ] P2 Finish the lark-based XPP parser path and configuration handling. refs: `clode/xpp_parser.py`
- [ ] P3 Export paths from internal RHS form back to Python/XPP. depends: stronger RHS IR. refs: `clode/xpp_parser.py`, `clode/function_converter.py`

## Performance, Randomness, And Scaling

- [ ] P1 Ensemble batching for device-capacity limits. depends: continuation semantics and integration/output separation. refs: `clode/trajectory.py`, `tmp/development_roadmap.md`
- [ ] P2 External RNG library audit (for example `Random123`) vs current RNG path. depends: reproducibility requirements. refs: `clode/kernels/clODE_random.cl`, `clode/kernels/features.cl`
- [ ] P2 Write-avoidance and memory-traffic audit for trajectory/event storage. depends: chunked/output redesign. refs: `clode/kernels/trajectory.cl`, `clode/kernels/odedriver.cl`

## Testing, Diagnostics, And Scope

- [ ] P1 Kernel-component tests with tiny synthetic models and kernels beyond end-to-end simulator tests. refs: `test/core_numerics/`, `clode/_opencl/source_builder.py`, `tmp/testing_audit.md`
- [ ] P1 Cache/build-key/invalidation coverage around program rebuilds and observer changes. refs: `clode/_opencl/program_cache.py`, `clode/_opencl/source_builder.py`, `clode/features.py`
- [ ] P2 Package/API hygiene: centralize defaults, typing cleanup, stepper-specific parameter subsets, trajectory output ergonomics. refs: `clode/solver.py`, `clode/features.py`, `clode/trajectory.py`
- [ ] P2 Publication-readiness and JOSS story: refresh the paper and repo landing page to the current PyOpenCL-only package, add a state-of-the-field comparison, and capture benchmark and impact evidence. refs: `paper/paper.md`, `README.md`, `tmp/joss_audit.md`
- [ ] P2 Scope audit: what belongs in clODE vs sibling/helper packages. refs: `tmp/package_state.md`, `tmp/development_roadmap.md`

## Inbox

- [ ] New rough idea with no clear home yet.
