# PyOpenCL Rollout Guardrails

## Purpose

This document records the rollout policy after feature parity landed.

It exists to answer three questions quickly:

- what validation bundles must stay green before changing rollout policy
- which runtime is currently considered stable on this workspace
- what must be true before switching the default backend

## Current Selector Policy

- the published package ships only the Python runtime path plus packaged OpenCL assets
- PyOpenCL is now the default backend for both tagged installs and source checkouts
- the C++ backend remains available only through explicit `_CLODE_BACKEND=cpp` selection when a locally built wrapper is present
- the C++ backend is now a comparison path rather than part of the default runtime story
- `pyopencl` is now part of the default package dependency set

## Stable Runtime On This Workspace

- the currently stable validation target is the NVIDIA runtime selected by `CLODE_TEST_PLATFORM_ID=0` and `CLODE_TEST_DEVICE_ID=0`
- the Intel CPU runtime remains unstable for the `localmax` feature rebuild path on this machine
- do not infer the stable runtime only from `clinfo -l`; on this workspace the CLI ordering does not match the stable backend-validation mapping that was verified in practice
- `clode.print_opencl()` now follows the active backend path: with `_CLODE_BACKEND=pyopencl` it reports the PyOpenCL-visible platform ordering rather than the legacy wrapper ordering

Recommended verification commands:

```bash
clinfo -l
```

and from Python:

```python
import clode
clode.print_opencl()
```

If platform ordering looks ambiguous, confirm the actual runtime selected by the backend before treating a device tuple as authoritative.

## Validation Tiers

### Authoritative backend-phase gate

- scope: fast backend migration gate
- current size: 55 tests
- current command:

```bash
CLODE_TEST_PLATFORM_ID=0 CLODE_TEST_DEVICE_ID=0 /home/fletcherpa/envs/clode/bin/python -m pytest test/core_numerics/test_transient.py test/core_numerics/test_trajectory.py test/core_numerics/test_features_basicall.py test/core_numerics/test_stochastic.py test/test_backend_contracts.py test/test_backend_rhs_source.py test/test_pyopencl_models.py test/test_pyopencl_source_builder.py test/test_pyopencl_runtime.py test/test_pyopencl_buffers.py test/test_pyopencl_structs.py test/test_pyopencl_transient_backend.py test/test_pyopencl_trajectory_backend.py test/test_pyopencl_feature_backend.py -q
```

### Extended reference bundle

- scope: authoritative gate plus higher-value reference coverage that is now worth carrying through rollout work
- current size: 77 tests
- current command:

```bash
CLODE_TEST_PLATFORM_ID=0 CLODE_TEST_DEVICE_ID=0 /home/fletcherpa/envs/clode/bin/python -m pytest test/core_numerics/test_transient.py test/core_numerics/test_trajectory.py test/core_numerics/test_features_basicall.py test/core_numerics/test_stochastic.py test/test_backend_contracts.py test/test_backend_rhs_source.py test/test_pyopencl_models.py test/test_pyopencl_source_builder.py test/test_pyopencl_runtime.py test/test_pyopencl_buffers.py test/test_pyopencl_structs.py test/test_pyopencl_transient_backend.py test/test_pyopencl_trajectory_backend.py test/test_pyopencl_feature_backend.py test/test_vdp.py test/test_features.py test/test_aux_values.py test/test_ornl_thompson_a1.py test/test_opencl_builtins.py test/test_runtime.py test/test_logger.py -q
```

This bundle passed on the stable NVIDIA runtime on this workspace after the PR 16 default-switch updates.

## Tests Still Excluded From PR12

- `test/test_solver.py` is not an active rollout target because its contents are mostly skipped placeholders
- the old empty `test/test_trajectory.py` file has been removed
- `test/test_observers.py` is manual/debug-oriented and currently skipped
- `test/test_clODE_utilities.py` is a placeholder note file, not an executable suite
- parser and converter tests remain useful, but they are not OpenCL-backend rollout blockers

## Default-Switch Prerequisites For PR16

Do not switch the default backend until all of the following are true:

1. the extended reference bundle is stable on the supported runtime set
2. the supported runtime policy is documented clearly enough that users can choose a working device tuple without guesswork
3. the PyOpenCL dependency and installation story are explicit enough for the default path
4. the C++ backend remains available as an explicit fallback during the transition
5. the backend selector and diagnostics make it obvious which runtime was actually chosen
6. the default package build no longer requires the C++ extension at import time

Status on this workspace after PR 16:

- satisfied: 1, 2, 3, 4, 5, and 6 for the current stable NVIDIA runtime
- PR 16 resolved the last selector-policy blocker: source checkouts no longer prefer the legacy backend merely because a local wrapper binary is present
- still unresolved for broader support: the Intel CPU runtime remains unstable for the `localmax` rebuild path

## Post-Switch Follow-Up Worth Keeping In Scope

- improve runtime-selection diagnostics so platform and device enumeration are easier to reconcile across `clinfo`, the legacy wrapper, and PyOpenCL
- separate supported-runtime policy from one-machine-specific observations once broader validation data exists
- revisit whether the Intel CPU runtime instability is a backend bug, a driver issue, or both before treating it as supported
