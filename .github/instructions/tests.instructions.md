---
applyTo: "test/**"
description: Use when editing or adding tests. Choose the smallest test layer that proves the intended contract.
---

# Test Placement Rules

- Use `test/core_numerics/` for exact-solution, convergence, and numerical-correctness evidence.
- Use `test/kernel_components/` for helper, build-key, storage, observer, and kernel-contract checks.
- Use the broader runtime and simulation tests for public API and behavior contracts, not for every internal plumbing assertion.
- Prefer direct component tests over broad end-to-end tests when they cover the same contract.
- Keep `tools/run_test_bundle.py` as the authoritative bundle map and avoid blurring gating correctness work with non-gating performance experiments.