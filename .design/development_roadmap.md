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
- The exact-solution and contract tests are much better than they were before the cleanup. The current package has a credible numerical regression spine.
- Split-window continuation correctness is now grounded by targeted regressions for `basicall` features and seeded stochastic Euler, and the live runtime has standardized on solver-owned absolute time plus persisted stochastic continuation details.
- The pure-Python packaging path and packaged kernel assets are now aligned with how a modern scientific Python package should ship.

### Current weaknesses

- Continuation correctness is much better, but the state model is still only partially explicit. Solver-owned time, attained per-item `tf`, continued `dt`, observer state, and RNG continuation details are real in the implementation, but they are not yet represented as one clear internal model.
- The current problem layer still stops short of an explicit IVP model. `ProblemInfo` and `RhsSource` exist, but the required default parameter and initial-state values still live in `Simulator` constructor arguments.
- Simulators still carry too much semantic weight. Ensemble generation, broadcasting, and some continuation semantics still live in `Simulator` helpers instead of first-class ensemble/state concepts.
- Each work item effectively owns `dt` and attained `tf`, but not an explicit `t0` or completion/error flags. That leaves friction around continuation helpers and windows where work items diverge.
- `SolverParams` still mixes integration policy with output-storage policy. That makes chunking, batching, and trajectory streaming harder than they need to be.
- Stepper, observer, solver-state, and problem concepts are still split across public enums/dataclasses, internal string registries, and kernel include trees rather than clearer Python-owned semantic definitions.
- Observer metadata is now Python-owned, but it is still expressed as a large conditional manifest rather than a cleaner observer-definition model.
- Optional event storage is still entangled with persistent observer state and build-time sizing.
- The runtime is intentionally single-device only. If multi-device execution ever becomes worthwhile, it will need a dedicated API and execution model rather than an extension of the current selectors.
- There is still no good path for stiff systems or implicit stepping, which limits the package for an important class of dynamical-systems problems.

## Prioritized Workstreams

### Priority 0: InitialValueProblem first, with built-in batch semantics at the simulator boundary

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
- add helper functions for grid, random, quasi-random, and repeat-style batch generation before deciding whether a separate `IVPEnsemble` class is warranted
- treat simulators as orchestration objects around IVP, runtime executor, and observer or trajectory policy
- keep explicit solver-state, observer-state, and stepper-definition cleanup as follow-on internal work unless a small adapter falls out naturally

Expected payoff:

- a clearer semantic unit for what one solve and one batch represent to the user without forcing a new container type prematurely
- simulator classes that read more cleanly as orchestration objects
- a cleaner base for later batching, continuation-policy helpers, and solver-state cleanup

API impact:

- can start internally; public API may only need helpers later

### Priority 1: Separate integration state from output and storage policy

This is the next structural bottleneck after continuation correctness.

Key tasks:

- decouple integration control from trajectory-storage capacity
- introduce chunked trajectory streaming or paged readback instead of treating `max_store` as one monolithic device allocation knob
- make room for batching large ensembles and long trajectories without forcing giant buffers
- separate persistent observer state from optional event-output capacity where possible

Why this matters:

- trajectory and feature storage policy currently leaks into core solver configuration
- large dynamical-systems scans will eventually need streaming and batching more than they need new public API surface

API impact:

- can start internally, but some cleanup of the public parameter story will probably come later

### Priority 1: Observer-definition and metadata cleanup

The current observer system works, but it is not yet shaped as a durable internal model.

Key tasks:

- replace the large observer-name conditionals with explicit Python-owned observer definitions
- model feature-name generation, persistent state layout, optional event layout, and warmup requirements in one place
- preserve the current kernels while making the metadata model easier to test and extend

Why this matters:

- observer infrastructure is a real differentiator for clODE
- future custom observers, event-storage refactors, and continuation fixes all benefit from a better observer-definition layer

Note on `odedriver.cl`:

- keep it as deferred design context for possible later unification, but do not use it as the starting point for the next PR sequence

### Priority 1: Explicit solver-owned state and stepper semantics

The lower-level execution model still needs a clearer internal contract, but it no longer has to define the active user-facing PR boundary.

Key tasks:

- introduce a per-work-item solver-state model with a clear home for `t0`, `tf`, `dt`, status flags, and continuation-specific RNG state
- use that model as the basis for continuation-policy helpers and later batching behavior where per-item windows can diverge
- introduce a Python-owned stepper-definition layer only after that state contract is clearer
- keep observer-definition cleanup adjacent but separate; `ObserverData` and optional event storage are their own modeling problem

Why after Priority 0:

- IVP-first batch cleanup clarifies what simulators orchestrate
- the solver-state layer can then be staged as an internal contract between `simulation` and `_opencl` instead of being entangled with user-facing IVP semantics

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

1. Introduce explicit IVP semantics and basic batch helpers while keeping simulator classes orchestration-focused.
2. Make lower-level solver state explicit, then use it as the contract for continuation helpers and later stepper cleanup.
3. Decouple integration state from output and storage policy so chunking and batching become feasible.
4. Clean up observer definitions and separate persistent observer state from optional event storage.
5. Tackle numerical kernel refinements such as fixed-step endpoint handling, `t0 + step * dt`, and interpolation improvements.
6. Do implicit-solver groundwork, then discuss the public API changes needed to expose it well.
