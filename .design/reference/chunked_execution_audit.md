# Chunked Execution Audit

Purpose: evaluate future time-chunking and ensemble-batching work, including the current capabilities, the likely benefits, and the internal boundaries that should stay clean now.
Read when: deciding how the landed solver-state boundary and the landed output/storage split should preserve room for chunked execution, progress reporting, long-run trajectory handling, or device-capacity batching.
Update when: continuation semantics, output-policy boundaries, batching strategy, or the recommended public API direction for chunked execution changes.

## Bottom line

- Treat time chunking and ensemble batching as two orthogonal execution-tiling policies layered on one semantic problem definition, not as new semantic owners.
- The current runtime already has the right low-level execution atom for a first implementation: repeated solves with explicit continuation of `x0`, `dt`, attained `tf`, RNG state, and observer state where applicable.
- The landed first-pass solver-state cleanup now makes that execution atom explicit and stable. Chunking still should not ship yet, but follow-on work should avoid re-coupling cache or buffer ownership in ways that would force another rewrite.
- The landed output/storage split is now the enabling layer that makes chunked trajectory reads and device-capacity batching straightforward rather than ad hoc.

## Current live capability

### What already works manually

- `Simulator.transient()` can already be repeated across windows while preserving device-side continuation state when `update_x0=True`.
- `FeatureSimulator.features()` can already be repeated while preserving observer state by default, so chunked feature accumulation is conceptually present today.
- `TrajectorySimulator.trajectory()` can already be repeated across windows, with each call returning only the current window so the caller can concatenate results on the host.
- Exact absolute-time continuation already has a verified manual pattern: run once, read `get_final_time()`, set the next `t_span` from the attained time, then run again.

### What is still awkward or missing

- There is no explicit chunking helper or execution policy; users have to hand-roll the loop.
- `shift_tspan()` advances the requested window, not the attained final time, so it is not the right primitive for exact chunked continuation.
- Public `SolverParams` still mixes integration controls with trajectory-storage controls (`max_store`, `nout`), but the internal and kernel-facing execution path now treats those as separate concerns.
- The ensemble still uses one shared requested `t_span`, so step-budget chunking for adaptive or otherwise diverged work items does not yet have a clear shared-window policy.
- Trajectory storage is still monolithic per call: `max_store` sizes one device allocation for the whole requested window instead of expressing a chunk or page policy.
- There is no orchestration layer for splitting a large IVP-owned ensemble into device-capacity-sized batches and then reassembling results.
- There are no progress hooks, cancellation hooks, or chunk-level output sinks yet.

## Why chunking and batching are worth doing

### Time chunking benefits

- Longer runs without requiring one giant trajectory or event buffer allocation.
- More interactive workflows: progress bars, intermediate inspection, early stopping, and notebook-friendly iteration.
- Cleaner handling of long fixed-step or stochastic runs where users want checkpoints or periodic summaries.
- Better control over host readback frequency instead of tying all output to one terminal solve call.
- A natural place to attach future progress reporting, checkpoints, or adaptive host-side policies.

### Ensemble batching benefits

- Run IVPs larger than device memory would allow in one monolithic launch.
- Keep trajectories or observer outputs bounded even when the full ensemble is very large.
- Make large parameter sweeps more resilient on smaller GPUs and more predictable in memory usage.
- Create a later path for balancing throughput against responsiveness without changing the semantic IVP model.

### Why these should be considered together

- They are different axes of the same execution-tiling problem: time chunks split one long solve along the time axis, while batching splits one large solve along the ensemble axis.
- A future orchestrator will likely need to compose them, for example by running multiple ensemble batches over multiple time chunks while preserving per-batch solver state.
- Both features depend on the same internal preconditions: explicit solver state, explicit output policy, and a clear distinction between semantic owners and transfer caches.

## Recommended mental model

Think in terms of five concerns that should stay separate:

- `InitialValueProblem`: the full semantic problem definition plus the full user-facing ensemble.
- solver state: per-work-item execution state for the active solve or batch, including requested window, attained `tf`, current `dt`, status, and RNG continuation.
- persistent observer state: feature-specific state that survives across chunks when the observer is not reinitialized.
- output policy: what to fetch, retain, page, or yield at each chunk boundary.
- execution tiling policy: how to split work over time chunks and ensemble batches.

Recommended rule:

- chunking and batching should orchestrate those owners; they should not become additional semantic owners themselves.

## Recommended API direction

### Internal direction

- Keep one shared internal chunk loop that reuses the existing repeated-solve execution atom instead of introducing a separate kernel mode first.
- Make chunk boundaries a Python-level orchestration concern first. Only move chunk loops into kernels later if profiling shows real benefit.
- Model three policy dimensions explicitly, even if they do not become public classes immediately:
  - continuation policy: requested-window vs attained-`tf`
  - chunk policy: stop after a time sub-interval or after a step budget
  - output policy: no fetch, final-state snapshot, trajectory page, or feature snapshot

### Public direction

- Keep the public API workflow-shaped rather than exposing a generic execution engine.
- If chunking becomes user-facing, prefer workflow-specific entry points backed by shared internals, such as chunked transient, trajectory, or feature helpers, rather than a new `ChunkedSimulator` type.
- Do not make users think about `_opencl` buffers, build keys, or cache objects.
- Keep continuation policy explicit. A later helper may default to attained-`tf` continuation, but it should not silently hide the distinction between requested and attained time.

## How time chunking likely differs by workflow

### Transient-only solves

- This is the easiest case.
- The main payload at each chunk boundary is updated solver state and optional final-state snapshots.
- A future progress API can likely hang off this path first.

### Feature accumulation

- This is already conceptually chunkable because observer state persists across repeated `features()` calls unless explicitly reinitialized.
- The key semantic question is whether chunk boundaries should expose cumulative observer outputs, per-chunk deltas, or both.
- Current evidence suggests the default should remain cumulative outputs unless a separate delta-oriented feature API is introduced later.

### Trajectory storage

- This is the main driver for separating integration state from output/storage policy.
- Each chunk naturally yields one trajectory page for the current requested sub-window.
- Host-side concatenation should remain explicit, including duplicated-boundary handling, unless a later sink or accumulator helper standardizes that policy.

## Relationship to ensemble batching

- The IVP should remain the semantic owner of the full ensemble and its shape metadata.
- Device-capacity batching should slice or view that full ensemble for execution; it should not force users to redefine the IVP as many smaller semantic problems.
- A later batching helper can live near the simulator or IVP-helper layer, but it should still preserve one full-ensemble semantic owner and one result-reassembly story.
- Chunking and batching should compose cleanly: the same run may need both a batch size and a chunk size.

Recommended rule:

- do not introduce a second top-level public ensemble owner just to make device batching possible.

## Implications for the current next PR

The next PR still should not implement chunking, but it should preserve the landed solver-state and output-policy contract.

### Decisions that help later chunking

- Give solver state one explicit internal home distinct from IVP-owned problem data, observer state, and fetched outputs.
- Include per-work-item `t0` in that solver-state model even if the first pass keeps it Python-owned rather than immediately threading it through every kernel.
- Keep attained `tf` and continued `dt` as first-class solver-state facts rather than incidental caches.
- Make runtime-only changes such as `t_span` updates and chunk boundaries avoid program rebuilds.
- Keep observer state persistent and adjacent to solver state, but not inside the solver-state object.
- Make fetched outputs derived caches, not the semantic owner of execution progress.
- Keep the solve loop legible in Python so later chunk helpers can layer on progress reporting and host-side policies without piercing `_opencl` internals.
- Do not force the first chunking-enabling PR to adopt a full matched device-side solver-state struct unless the narrower state cleanup proves insufficient.

### Decisions that would make later chunking harder

- Treating `max_store` or `nout` as part of core solver semantics instead of output policy.
- Making `shift_tspan()` or another convenience helper the hidden owner of continuation policy.
- Letting simulator-side `_device_*` caches continue to stand in for explicit solver state.
- Tying observer persistence too tightly to feature-output allocation details.

## Recommended follow-on order

1. Keep the landed solver-state and output-policy cleanup as the contract.
2. Clean up observer definitions so feature chunking does not have to fight ambiguous observer state/layout boundaries.
3. Add one internal chunk-orchestration path for repeated solves and progress hooks.
4. Add trajectory paging or chunked readback helpers on top of that output-policy layer.
5. Add device-capacity ensemble batching that composes with the same chunk loop.
6. Add public ergonomic helpers once the internal policies stop moving.

## Test implications for later work

- Split-window chunked transient runs should match one long run for deterministic and seeded stochastic cases.
- Split-window chunked feature runs should match one long run when continuation uses attained `tf` and observer state persists.
- Chunked trajectory pages should stitch back together without semantic drift beyond the documented duplicated-boundary rule.
- Results should be invariant to ensemble batch size.
- Chunk or batch boundaries should not trigger unnecessary rebuilds or break program-cache reuse.

## Open questions to keep explicit

- Should step-budget chunking become public before the package has a clearer per-work-item status model?
- Should chunked feature APIs expose cumulative observer state only, or also per-chunk deltas?
- Should public batching appear first as a simulator option, an IVP helper, or a standalone orchestration helper?
- How much progress and checkpointing support belongs in clODE itself versus notebook or application code layered on top?
