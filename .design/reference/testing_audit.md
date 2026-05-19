# Testing Audit

## Bottom line

The current bundle and marker taxonomy is a good base. The new `test/kernel_components/` layer now provides the missing middle layer between full end-to-end simulator tests and the small host-side `_opencl` support tests, but its current coverage is intentionally narrow and still needs expansion.

## What is already working well

- `tools/run_test_bundle.py` already gives the suite a useful domain taxonomy: `smoke`, `frontend`, `runtime_api`, `numerics`, `kernel_components`, and `opencl_internal`.
- `test/conftest.py` centrally injects markers by file and path, which keeps selection logic out of individual test modules.
- `test/core_numerics/` already has the right shape for authoritative numerical regressions: small models, exact references, and reusable helpers.
- `test/kernel_components/` now directly covers helper-kernel behavior and a basic-observer contract without going through the full simulator stack.
- CI already separates a broad cross-platform smoke matrix from a narrower OpenCL-backed release gate.

## Current weaknesses

- The physical layout is still transitional, so the intent of each file is not always obvious from its location.
- Too much OpenCL-specific correctness signal still arrives through full simulator paths instead of smaller component tests; the new component layer is still only a first slice.
- There is no clearly separated performance layer yet.
- Device and platform coverage is still mostly controlled by environment selection instead of a richer parametrization story.

## Recommended test taxonomy

Keep the current bundles, but gradually organize tests into these layers:

- `smoke`: import, packaging, minimal frontend checks, and anything that must stay driver-independent.
- `frontend`: parser, function-converter, builtins, and user-facing ingestion paths.
- `runtime_api`: public runtime selection, simulator contracts, defaults, logging, and backend-facing public behavior.
- `core_numerics`: authoritative correctness tests against exact solutions or fixed expectations.
- `kernel_components` or `opencl_components`: the middle layer for build options, struct layout, buffer semantics, observer storage behavior, and compile-time feature toggles.
- `opencl_internal`: host-side support-layer tests that are genuinely internal to the OpenCL implementation.
- `performance`: explicitly non-gating throughput, allocation, and compile-latency checks.

The key point is not a perfect directory tree. The key point is separating correctness, API, component, and performance intent more clearly.

## Correctness vs performance vs API vs internal

Use these rules:

- Correctness tests answer whether the numerical or semantic result is right.
- API tests answer whether users can call the package the documented way and get the documented behavior.
- Internal and component tests answer whether a specific buffer, build, struct, or kernel contract still holds.
- Performance tests answer whether the package remains fast enough on a known device class, but they should be signal-only unless the project later commits to stable benchmark hardware.

## Custom kernel strategy

clODE should add tiny synthetic kernels and tiny synthetic models on purpose instead of relying only on larger biological or workflow-driven models.

Good uses for custom kernels:

- Struct layout and buffer-shape checks.
- Build-key invalidation and observer rebuild checks.
- Observer initialization and finalization paths.
- Event-storage sizing and `N_STORE_EVENTS` behavior.
- Utility-kernel behavior that does not need a full simulator run to validate.

Good uses for synthetic models:

- Exact-solution transient and trajectory regressions.
- Minimal auxiliary-variable coverage.
- Reproducibility and continuation semantics.
- Boundary-case models that force a specific observer or stepper path.

Rule of thumb:

- Use analytic or reference models when testing numerical correctness.
- Use tiny synthetic kernels when testing plumbing, storage, or build behavior.

## CI guidance

- Keep the current `release = frontend + runtime_api + numerics` gate narrow and meaningful.
- Do not put `opencl_internal` or `performance` into the default release gate.
- Add a new component-oriented bundle first for manual or Linux OpenCL validation, then decide later whether parts of it belong in the gate.
- Consider PyOpenCL's own pytest parametrization helpers for signal-only multi-device coverage.

## Concrete next steps

1. Expand `test/kernel_components/` beyond helper math and the `basic` observer path.
2. Move future build-key, observer-storage, and source-assembly checks there instead of burying them inside larger end-to-end tests.
3. Add a non-gating performance bundle with explicit hardware assumptions.
4. Keep `tools/run_test_bundle.py` authoritative even if the physical file layout changes slowly.

## Recommendation

Do not throw away the current bundle taxonomy. Build on it by adding a component-test layer that uses minimal kernels and minimal models to test OpenCL-specific contracts directly.
