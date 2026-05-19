# clODE Development Roadmap

Purpose: medium-lived rationale and priority ordering for follow-on work.
Read when: you need tradeoffs, architectural sequencing, or the why behind a workstream.
Update when: priorities shift materially or the architectural rationale changes.

Use `.design/ideas.md` as the short living board. This file is the longer rationale and prioritization note.

## Evaluation

### Current strengths

- The numerical core is still strong: one work-item per trajectory, compile-time specialization, and a clean separation between stepper kernels and observer kernels.
- The PyOpenCL migration succeeded in the important place: runtime ownership, source assembly, build keys, and program caching now live in Python rather than in a compiled wrapper.
- `FeatureSimulator` remains a distinctive strength. Stateful observers let the package extract oscillation and event statistics without forcing full trajectory storage.
- Observer definitions, resolved specs, and observer-state naming now line up across Python, OpenCL metadata, and the active kernels, which makes feature builds and invalidation rules much easier to reason about.
- The exact-solution and contract tests are much better than they were before the cleanup. The current package has a credible numerical regression spine.
- Split-window continuation correctness is now grounded by targeted regressions for `basicall` features and seeded stochastic Euler, and the live runtime has standardized on solver-owned absolute time plus persisted stochastic continuation details.
- The pure-Python packaging path and packaged kernel assets are now aligned with how a modern scientific Python package should ship.

### Current weaknesses

- Continuation correctness is much better, and the first-pass Python-owned solver-state model, the integration/output split, the observer-definition cleanup, the execution-setting cleanup, and the stepper-definition cleanup are now live all the way into the OpenCL layer. The next structural bottleneck is continuation ergonomics: the solver owns absolute time, but callers still have to express requested-window versus attained-`tf` policy by hand.
- A first-pass IVP model now exists, and the curated public problem layer now centers `InitialValueProblem`; lower-level support types and source-preparation helpers are internal support concepts rather than promoted surface area.
- Simulators still carry some semantic weight through compatibility delegates and cached mirrored problem data even though the core defaults and batch semantics now live on the IVP.
- The next internal friction point is clearer now: the internal execution model is much more coherent, but repeated-solve continuation still feels lower-level than it should because callers must manage attained-`tf` versus requested-window policy explicitly.
- Each work item effectively owns `dt` and attained `tf`, but not an explicit `t0` or completion/error flags. That leaves friction around continuation helpers and windows where work items diverge.
- `SolverParams` still mixes integration policy with output-storage policy. That makes chunking, batching, and trajectory streaming harder than they need to be.
- Stepper, solver-state, and problem concepts are still split across public enums/dataclasses, internal string registries, and kernel include trees rather than clearer Python-owned semantic definitions.
- The runtime is intentionally single-device only. If multi-device execution ever becomes worthwhile, it will need a dedicated API and execution model rather than an extension of the current selectors.
- There is still no good path for stiff systems or implicit stepping, which limits the package for an important class of dynamical-systems problems.

## Prioritized Workstreams

### Priority 0: InitialValueProblem first, with built-in batch semantics at the simulator boundary

Status on the current branch:

- the first semantic IVP pass has landed
- simulators can now consume `ivp=` directly while preserving the compatibility constructor path
- Python-authored IVPs can act as SciPy-style `(t, y, *args) -> dydt` callables
- the remaining work in this priority is now mostly future batch-helper follow-through rather than more boundary expansion in the same PR

Continuation correctness is no longer the headline problem. The next durable improvement is to make the user-facing semantic unit explicit where the current API already points: one IVP, which can also cover the usual size-`(1,)` case plus basic batched inputs.

Why now:

- the public API already treats default parameter and initial-state values as mandatory, which makes `InitialValueProblem` a more natural first semantic object than a separate `ProblemDefinition`
- `set_ensemble()` and `set_repeat_ensemble()` already expose a distinct batch-shaping concern, but that concern still lives as simulator-side array shaping and broadcasting logic
- simulator classes currently mix IVP/default construction, ensemble shaping, and solve orchestration in one place
- an IVP that owns basic batch semantics and remembered shape may be enough for the first pass, which is likely cheaper for users than introducing a dedicated public ensemble class too early
- clarifying the simulator boundary first should make later solver-state, continuation-policy, and stepper work easier to stage without forcing them into the same PR

Key tasks:

- introduce an `InitialValueProblem` model that owns RHS semantics plus default parameter and initial-state values
- keep any lower-level shared definition metadata derived or internal unless a separate layer proves necessary later
- let the IVP own the size-`(1,)` case plus basic batched parameter or initial-state semantics, with remembered shape metadata where that improves result reshaping workflows
- ensure an intuitive user experience for IVP specification. Evaluate callable IVPs for SciPy `solve_ivp`, but bias toward a small Python-backed adapter on the IVP rather than widening this PR into generic OpenCL/XPP round-tripping or converter redesign
- add helper functions for grid, random, quasi-random, and repeat-style batch generation before deciding whether a separate `IVPEnsemble` class is warranted
- treat simulators as orchestration objects around IVP, runtime executor, and observer or trajectory policy
- keep explicit solver-state, observer-state, and stepper-definition cleanup as follow-on internal work unless a small adapter falls out naturally

Expected payoff:

- a clearer semantic unit for what one solve and one batch represent to the user without forcing a new container type prematurely
- simulator classes that read more cleanly as orchestration objects
- a cleaner base for later batching, continuation-policy helpers, and solver-state cleanup

API impact:

- can start internally; public API may only need helpers later. If callable IVPs land here, keep them conditional on retaining a Python RHS callable on the IVP and defer generic cross-source interop until stronger RHS IR work.

### Landed Priority 0: Separate integration state from output and storage policy

Status on the current branch:

- the split is now explicit in simulator settings views, executor invalidation rules, OpenCL buffers, matched struct packing, and kernel entrypoints
- common runtime state now carries integration settings only
- trajectory buffers own trajectory output policy, and runtime observer settings no longer include event timestamp capacity
- trajectory and feature executors now guard against stale direct downloads when output-policy changes invalidate prior results

What this unlocked:

- chunked trajectory paging and device-capacity batching no longer have to fight a mixed solver/output ABI first
- observer work can now focus on observer definitions and state/layout boundaries rather than reopening solver/output ownership
- broader public config redesign can stay deferred until the observer and stepper semantics stop moving

### Landed Priority 0: Single source of truth for execution-setting defaults and compatibility resolution

Status on the current branch:

- canonical default integration and trajectory-output settings now live in one place under `clode/simulation/params.py`
- `SolverParams` defaults and simulator constructor defaults no longer drift
- transient, trajectory, and feature simulators now resolve scalar solver arguments and prebuilt `SolverParams` bundles through the same helper
- simulators keep internal copies of caller-provided `SolverParams` bundles instead of aliasing them

What this unlocked:

- stepper work can now build on stable execution-setting semantics instead of another layer of duplicated defaults
- later public config redesign can stay focused on API shape rather than cleanup of internal default-resolution drift

### Landed Priority 0: Observer-definition and observer-state cleanup

Status on the current branch:

- explicit Python-owned observer definitions and resolved specs now drive feature names, build defines, layout boundaries, and feature-executor invalidation
- persistent observer state is now distinct from runtime settings and from optional event-output capacity
- the active kernels, matched struct lookups, and synthetic host-side names now use observer-state naming consistently

What this unlocked:

- later stepper and config work no longer need to reopen observer layout or build semantics
- future custom observers and any later observer-specific public API discussion can build on a clearer internal model first

Note on `odedriver.cl`:

- keep it as deferred design context for possible later unification, but do not use it as the starting point for the next PR sequence

### Priority 0: Explicit solver-owned state and cache ownership cleanup

Status on the current branch:

- the first pass has landed
- `clode/simulation/_state.py` now owns Python-side solver state and fetched-output caches
- IVP owns next-solve problem data while `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners
- direct invalidation and continuation regressions now cover the new boundary

What remains for later:

- no device-side per-work-item `t0` or richer completion/error status model yet
- exact absolute-time continuation still needs caller-managed `t_span` or a later public continuation helper
- the next leverage point is a public continuation-policy helper rather than reopening the ownership or output-policy work itself

### Landed Priority 0: Python-owned stepper semantics after execution-setting cleanup

Status on the current branch:

- built-in stepper traits and OpenCL build mapping now resolve through `clode/simulation/_stepper_definitions.py`
- `_opencl/registry.py`, `_opencl/source_builder.py`, and executor construction now consume the stepper-definition catalog instead of raw string tables
- build keys now carry the resolved stepper define so program-cache identity tracks the effective kernel specialization more directly

What this unlocked:

- continuation helpers and later implicit-stepper groundwork can build on a clearer internal semantic layer first
- source assembly and runtime validation now have one coherent home for stepper semantics instead of split registry logic

### Priority 0: Public continuation-policy helper/API after stepper cleanup

The next active cleanup should turn the solver-owned time model into a more explicit, ergonomic public policy.

Why this should be next:

- the internal solver-state and stepper-definition boundaries are now stable enough that continuation ergonomics no longer need to sit on shifting internals
- the current need to hand-roll `set_tspan(get_final_time())` is both easy to misuse and harder to explain than it should be
- this is a high-value developer and user clarity improvement without forcing a broader public config redesign yet

The lower-level execution model still needs a clearer internal contract, but it no longer has to define the active user-facing PR boundary.

Key tasks:

- add one explicit continuation-policy helper or small API that distinguishes requested-window continuation from attained-`tf` continuation
- keep the helper aligned with the current solver-owned time semantics and current shared-`tspan` runtime limitations
- make the common continuation path more discoverable without hiding the remaining per-work-item `t0` limitations

Why after Priority 0:

- IVP-first batch cleanup clarifies what simulators orchestrate
- the solver-state layer and observer-definition layer are now clearer internal contracts between `simulation` and `_opencl`
- execution-setting cleanup and stepper-definition cleanup have now landed, so continuation ergonomics can build on stable internal semantics instead of another moving target

### Priority 2: Numerical kernel refinements

These items are important, but they should follow the state-model work.

Key tasks:

- address fixed-step endpoint drift without violating fixed-step semantics
- consider more accurate time accumulation helpers such as `TwoSum`, `t0 + step * dt`, or related compensated/structured-time ideas already noted in the kernels
- improve interpolation and dense-output groundwork for trajectory and observer use cases

Why later:

- these are meaningful numerical improvements, but continuation mismatches are a more immediate correctness issue

### Priority 2: Solver-extension groundwork for stiff problems

Implicit methods are a natural long-term direction, but they should not be the immediate next implementation item.

Near-term groundwork:

- clarify how stepper definitions should carry method metadata beyond the current explicit/fixed-adaptive set
- decide what internal Jacobian and residual model is acceptable for OpenCL execution
- separate solver-state questions from output and observer concerns before adding nonlinear solver machinery

Why not first:

- the current package still has unresolved state and continuation semantics
- adding implicit methods now would build on abstractions that are not yet stable enough

API impact:

- likely blocked from full public exposure until collaborators are ready for API discussion

### Priority 3: Runtime-surface and package cleanup

This is worthwhile, but it is lower leverage than the items above.

Key tasks:

- decide whether the package should eventually keep only `clode.__init__` as the stable top-level barrel and retire the remaining flat compatibility modules through a normal deprecation cycle
- keep future runtime/package cleanup scoped to real API clarity rather than historical naming debt

Why later:

- these are architecture-hygiene improvements, not the main scientific bottlenecks right now

Packaging note:

- the flat compatibility barrels are not needed for package discovery or wheel correctness
- their only real value is import-path continuity and module-path stability for downstream users during the transition

## Longer-Term Items That Likely Need API Discussion

- true multi-device execution via a dedicated API rather than overloading the current single-device selectors
- a public observer-definition or custom-observer interface
- explicit public solver-state or observer-state objects
- public exposure of implicit or linearly implicit steppers
- a cleaner separation between integration parameters and output-policy parameters in the public constructors

## sdist Note

- Omitting docs from the sdist is reasonable if the goal is a smaller, cleaner source artifact.
- Keeping tests in the sdist is still the better default unless artifact size becomes a real problem, because downstream packagers and maintainers benefit from having the verification surface available.

## Recommended Order Of Attack

1. Add a public continuation-policy helper on top of the landed solver-state and stepper-definition boundaries.
2. Tackle numerical kernel refinements such as fixed-step endpoint handling, `t0 + step * dt`, and interpolation improvements.
3. Do implicit-solver groundwork.
4. Revisit broader public config/API redesign only after continuation ergonomics and stepper semantics are stable.

Public continuation helpers, richer IVP batch helpers, and broader solver interop should follow once those internal boundaries are stable enough that they are unlikely to be redesigned immediately afterward.
