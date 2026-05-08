# clODE Development Roadmap

Use `tmp/ideas.md` as the short living board. This file is the longer rationale and prioritization note.

## Evaluation

### Current strengths

- The numerical core is still strong: one work-item per trajectory, compile-time specialization, and a clean separation between stepper kernels and observer kernels.
- The PyOpenCL migration succeeded in the important place: runtime ownership, source assembly, build keys, and program caching now live in Python rather than in a compiled wrapper.
- `FeatureSimulator` remains a distinctive strength. Stateful observers let the package extract oscillation and event statistics without forcing full trajectory storage.
- The exact-solution and contract tests are much better than they were before the cleanup. The current package has a credible numerical regression spine.
- The pure-Python packaging path and packaged kernel assets are now aligned with how a modern scientific Python package should ship.

### Current weaknesses

- Continuation semantics are still only partially modeled. Solver state, observer state, requested time windows, attained final times, and cached outputs are related, but they are not yet represented as a clear internal state model.
- Two important continuation issues still reproduce on the current PyOpenCL path:
  - `basicall` split-window feature continuation does not match a single long run.
  - seeded stochastic Euler split-window continuation does not match a single seeded run.
- `SolverParams` still mixes integration policy with output-storage policy. That makes chunking, batching, and trajectory streaming harder than they need to be.
- Observer metadata is now Python-owned, but it is still expressed as a large conditional manifest rather than a cleaner observer-definition model.
- Optional event storage is still entangled with persistent observer state and build-time sizing.
- The public runtime surface still advertises `device_ids`, even though the implementation is intentionally single-device only.
- There is still no good path for stiff systems or implicit stepping, which limits the package for an important class of dynamical-systems problems.

## Prioritized Workstreams

### Priority 0: Continuation correctness and explicit state semantics

This is the highest-value next area.

Why now:

- it is already failing in live reproduced cases
- it sits underneath batching, chunking, non-autonomous correctness, and future solver expansion
- it can be addressed without changing the public API

Key tasks:

- add permanent regressions for split-window feature continuation and seeded stochastic continuation
- make the code distinguish more clearly between requested `t_span` and attained final time
- model the continued solver state explicitly enough to reason about `x0`, `dt`, RNG state, and cached outputs separately
- model observer continuation state explicitly enough to reason about persistent `ObserverData`, initialization state, and feature readback separately

Expected payoff:

- safer repeated-call semantics
- better non-autonomous behavior
- a cleaner base for later batching, streaming, and implicit methods

API impact:

- none required if kept internal

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

### Priority 2: Numerical kernel refinements

These items are important, but they should follow the state-model work.

Key tasks:

- address fixed-step endpoint drift without violating fixed-step semantics
- consider more accurate time accumulation helpers such as `TwoSum` or related compensated summation ideas already noted in the kernels
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

- decide whether `device_ids` should eventually be removed or made real
- decide whether the package should eventually keep only `clode.__init__` as the stable top-level barrel and retire the remaining flat compatibility modules through a normal deprecation cycle
- keep future runtime/package cleanup scoped to real API clarity rather than historical naming debt

Why later:

- these are architecture-hygiene improvements, not the main scientific bottlenecks right now

Packaging note:

- the flat compatibility barrels are not needed for package discovery or wheel correctness
- their only real value is import-path continuity and module-path stability for downstream users during the transition

## Longer-Term Items That Likely Need API Discussion

- true multi-device execution rather than nominal multi-device arguments
- a public observer-definition or custom-observer interface
- explicit public solver-state or observer-state objects
- public exposure of implicit or linearly implicit steppers
- a cleaner separation between integration parameters and output-policy parameters in the public constructors

## sdist Note

- Omitting docs from the sdist is reasonable if the goal is a smaller, cleaner source artifact.
- Keeping tests in the sdist is still the better default unless artifact size becomes a real problem, because downstream packagers and maintainers benefit from having the verification surface available.

## Recommended Order Of Attack

1. Fix split-window continuation correctness and clarify internal continuation state.
2. Decouple integration state from output and storage policy so chunking and batching become feasible.
3. Clean up observer definitions and separate persistent observer state from optional event storage.
4. Tackle numerical kernel refinements such as fixed-step endpoint handling and interpolation improvements.
5. Do implicit-solver groundwork, then discuss the public API changes needed to expose it well.
