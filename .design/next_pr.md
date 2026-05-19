# Next PR

Purpose: the single active implementation target.
Read when: you need the current highest-priority PR scope, rationale, or acceptance criteria.
Update when: the active target changes, the scope narrows or broadens, or the acceptance criteria change.

## Title

Separate integration state from output and storage policy

## Assumed Repo State

- Canonical public homes are now `clode.problem`, `clode.observers`, `clode.simulation`, and `clode.runtime`.
- `clode._opencl` is the canonical internal execution layer.
- A first-pass internal solver-state boundary now exists: IVP owns next-solve problem data, `clode/simulation/_state.py` owns Python-side solver state and fetched-output caches, and `_opencl/executors.py` treats host mirrors as transfer caches rather than semantic owners.
- Split-window continuation regressions and the new ownership/invalidation regressions are green on the live PyOpenCL path.
- Root flat modules such as `clode.solver` and `clode.features` are compatibility barrels only; new work should target canonical packages unless the task is explicitly about compatibility cleanup.

## Why this should be next

The solver-state cleanup landed the boundary we needed, but it also made the next bottleneck more obvious: integration policy and output or storage policy are still coupled in too many places.

Today `SolverParams` still mixes integration controls with trajectory storage controls such as `max_store` and `nout`, and trajectory or feature buffer sizing still reads like part of solver state rather than an explicit output policy. That coupling is now the main reason chunking, batching, observer-state cleanup, and more ergonomic continuation helpers still feel awkward.

This is the highest-value next PR because it builds directly on the new solver-state contract without reopening that ownership work, and it is the cleanest internal prerequisite for chunked trajectory paging, device-capacity ensemble batching, and clearer observer-state versus event-storage semantics.

## Scope

- introduce an explicit internal boundary between integration state and output or storage policy
- keep solver state focused on requested window, continued time, attained `tf`, current `dt`, and continuation-related sync facts
- stop treating `max_store`, `nout`, and related feature or event output capacity as core solver-state facts
- make trajectory and feature/event buffer lifecycle follow explicit output-policy changes rather than broad solver invalidation
- preserve the current solver-state, IVP, and persistent observer-state ownership split while making output caches and output allocation rules thinner and more legible
- keep public API changes minimal; prefer internal cleanup first and compatibility-preserving adapters if a small public helper falls out naturally

## Design Constraints

- keep `InitialValueProblem` as the semantic owner of next-solve problem data only
- keep the landed `SolverState` boundary intact; do not move output policy back into solver state
- treat persistent observer state as adjacent runtime state, but keep optional event-output capacity outside the core solver-state model
- keep fetched outputs and transfer caches as derived data, not semantic owners
- separate build-affecting observer layout changes from runtime-only output-policy changes when possible
- preserve the current kernel-facing buffer ABI unless a narrow change materially clarifies the integration/output split
- keep chunking and ensemble batching as follow-on orchestration policies layered on these owners, not as new semantic owners

## Non-goals

- no new public IVP expansion beyond incidental polish
- no chunked trajectory streaming or ensemble batching yet
- no public continuation-policy helper in the same PR unless a tiny boundary-preserving helper falls out naturally
- no observer-definition redesign beyond the internal boundary cleanup needed to separate persistent observer state from optional output capacity later
- no device-side per-work-item `t0` or full matched solver-state ABI redesign in the same PR
- no multi-device work

## Suggested Implementation Slices

1. Introduce explicit internal output/storage policy types or helpers and route trajectory capacity plus feature/event allocation through them without changing numerical behavior.
2. Refactor simulator and executor invalidation so integration changes clear results without looking like output-policy mutations, while output-policy changes only rebuild or reallocate what they actually own.
3. Split the current `SolverParams` semantics internally so integration controls and output/storage controls stop moving together even if the public constructor stays compatible for now.
4. Verify the same separation across transient, trajectory, and feature paths, especially where persistent observer state and optional event storage interact.
5. Add focused regressions for invalidation, buffer reuse, and continued numerical correctness before any follow-on chunking or batching helper work.

## Code-Facing Checklist

- `clode/simulation/params.py`: clarify or stage the internal split between integration controls and output/storage controls while preserving the current public surface as needed.
- `clode/simulation/base.py`: keep solver-state invalidation scoped to solver concerns and avoid reintroducing output-policy ownership there.
- `clode/simulation/trajectory.py`: make trajectory fetch and reshape logic depend on explicit output-policy and trajectory-cache boundaries rather than on broad solver resets.
- `clode/simulation/features.py`: keep persistent observer state separate from feature-output invalidation and optional event-capacity changes.
- `clode/_opencl/buffers.py`: make the logical split between common solver buffers, trajectory storage, and feature/observer output allocation easier to follow.
- `clode/_opencl/executors.py`: keep transfer-cache invalidation grouped by concern and make trajectory or feature buffer reallocation follow output-policy changes rather than unrelated solver resets.
- `clode/_opencl/observer_metadata.py` and related feature-buffer code: only touch where needed to separate persistent observer state from optional output capacity.
- `test/test_simulation_contracts.py`: add or adjust ownership and invalidation coverage where solver state and output policy now diverge.
- `test/core_numerics/test_trajectory.py`: keep trajectory storage behavior and `nout`/`max_store` semantics under direct regression coverage.
- `test/core_numerics/test_features_basicall.py` and `test/core_numerics/test_stochastic.py`: keep continuation behavior green while the output-policy cleanup lands.

## Acceptance Criteria

- there is one explicit internal home for integration state and one explicit internal home for output/storage policy; they no longer read like one mixed concern
- `max_store`, `nout`, and feature/event output capacity no longer act like core solver-state facts in the internal code paths
- trajectory and feature/event buffer allocation follows output-policy changes rather than unrelated solver resets
- the landed solver-state boundary remains intact: IVP-owned problem data, solver-owned execution state, persistent observer state, and fetched outputs stay clearly separated
- split-window continuation regressions remain green for transient, `basicall`, and seeded stochastic cases
- the current public constructor and compatibility paths remain green unless a deliberate, documented change is made

## Follow-on If This Lands Cleanly

The next few high-value PRs should be observer-definition cleanup around persistent state versus optional event storage, then a cleaner Python-owned stepper-definition model, then chunked trajectory or batching helpers built on the separated output-policy layer. Public continuation helpers and richer IVP batch helpers should follow once those boundaries stop moving.
