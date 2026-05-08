# Next PR

## Title

Split-window continuation correctness

## Why this should be next

Two archived continuation issues still reproduce on the current PyOpenCL implementation on the stable workspace runtime (`CLODE_TEST_PLATFORM_ID=0`, `CLODE_TEST_DEVICE_ID=0`):

- `basicall` feature continuation does not match a single-window reference run. The reproduced max feature difference is about `0.26948` on `stable_linear_aux.cl` with RK4 and `dt=0.125`.
- seeded stochastic Euler split-window continuation does not match a single seeded run. The reproduced max state difference is about `0.84252` on the Ornstein-Uhlenbeck ensemble with `dt=0.125`.

These are not cosmetic issues. They sit directly on top of solver continuation state, observer continuation state, RNG semantics, and time-window semantics. If they remain unresolved, later work on chunking, batching, non-autonomous systems, and implicit solvers will be built on an unstable base.

## Scope

- add permanent OpenCL-backed regressions for:
  - split-window `basicall` feature continuation
  - split-window seeded stochastic Euler continuation
- trace and fix the stochastic boundary-state semantics so split windows consume the same effective noise sequence as a single uninterrupted run
- trace and fix the `basicall` continuation path so accumulated observer statistics are associative across windows
- make the internal continuation model explicit enough that the code distinguishes:
  - requested `t_span`
  - attained final time
  - continued solver state (`x0`, `dt`, RNG)
  - continued observer state
- keep the public API unchanged

## Non-goals

- no chunked trajectory streaming
- no implicit solver addition
- no multi-device work
- no public API redesign
- no broad package or module renaming

## Acceptance Criteria

- the new regressions fail on the current code and pass after the fix
- a split run that feeds the first segment's attained `tf` into the second segment matches the single-run reference within agreed tolerances for both reproduced cases
- existing release and numerics bundles remain green on the stable runtime
- the PR leaves behind clearer internal structure or comments around continuation state instead of only patching two isolated symptoms

## Follow-on If This Lands Cleanly

The next high-value step should be to separate integration state from output and storage policy (`max_store`, event storage, chunking, batching), because that refactor becomes much safer once continuation semantics are trustworthy.
