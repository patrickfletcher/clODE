# PyOpenCL Rollout Guardrails

## Purpose

This document records the rollout policy after feature parity landed.

It exists to answer three questions quickly:

- what validation bundles must stay green before changing rollout policy
- which runtime is currently considered stable on this workspace
- what must be true before switching the default backend

## Current Selector Policy

- the public default backend remains the current C++ path
- the PyOpenCL backend remains behind the internal `_CLODE_BACKEND=pyopencl` selector
- the C++ backend remains the required fallback while PyOpenCL rollout is still runtime-qualified
- the optional dependency path for the PyOpenCL backend is now exposed as `clode[pyopencl]`
- packaging audit note: the current default distribution is still a Bazel-backed C++ wheel, so a backend default switch alone does not finish the migration

## Stable Runtime On This Workspace

- the currently stable validation target is the NVIDIA runtime selected by `CLODE_TEST_PLATFORM_ID=0` and `CLODE_TEST_DEVICE_ID=0`
- the Intel CPU runtime remains unstable for the `localmax` feature rebuild path on this machine
- do not infer the stable runtime only from `clinfo -l`; on this workspace the CLI ordering does not match the stable backend-validation mapping that was verified in practice

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
- current size: 54 tests
- current command:

```bash
CLODE_TEST_PLATFORM_ID=0 CLODE_TEST_DEVICE_ID=0 /home/fletcherpa/envs/clode/bin/python -m pytest test/core_numerics/test_transient.py test/core_numerics/test_trajectory.py test/core_numerics/test_features_basicall.py test/core_numerics/test_stochastic.py test/test_backend_contracts.py test/test_backend_rhs_source.py test/test_pyopencl_models.py test/test_pyopencl_source_builder.py test/test_pyopencl_runtime.py test/test_pyopencl_buffers.py test/test_pyopencl_structs.py test/test_pyopencl_transient_backend.py test/test_pyopencl_trajectory_backend.py test/test_pyopencl_feature_backend.py -q
```

### Extended reference bundle

- scope: authoritative gate plus higher-value reference coverage that is now worth carrying through rollout work
- current size: 74 tests
- current command:

```bash
CLODE_TEST_PLATFORM_ID=0 CLODE_TEST_DEVICE_ID=0 _CLODE_BACKEND=pyopencl /home/fletcherpa/envs/clode/bin/python -m pytest test/core_numerics/test_transient.py test/core_numerics/test_trajectory.py test/core_numerics/test_features_basicall.py test/core_numerics/test_stochastic.py test/test_backend_contracts.py test/test_backend_rhs_source.py test/test_pyopencl_models.py test/test_pyopencl_source_builder.py test/test_pyopencl_runtime.py test/test_pyopencl_buffers.py test/test_pyopencl_structs.py test/test_pyopencl_transient_backend.py test/test_pyopencl_trajectory_backend.py test/test_pyopencl_feature_backend.py test/test_vdp.py test/test_features.py test/test_aux_values.py test/test_ornl_thompson_a1.py test/test_opencl_builtins.py test/test_runtime.py -q
```

This bundle passed on the stable NVIDIA runtime on this workspace, including after the PR 13 transition updates.

## Tests Still Excluded From PR12

- `test/test_solver.py` is not an active rollout target because its contents are mostly skipped placeholders
- `test/test_trajectory.py` is empty
- `test/test_observers.py` is manual/debug-oriented and currently skipped
- `test/test_clODE_utilities.py` is a placeholder note file, not an executable suite
- parser and converter tests remain useful, but they are not OpenCL-backend rollout blockers

## Default-Switch Prerequisites For PR16

Do not switch the default backend until all of the following are true:

1. the 74-test extended reference bundle is stable on the supported runtime set
2. the supported runtime policy is documented clearly enough that users can choose a working device tuple without guesswork
3. the PyOpenCL dependency and installation story are explicit enough for the default path
4. the C++ backend remains available as an explicit fallback during the transition
5. the backend selector and diagnostics make it obvious which runtime was actually chosen
6. the default package build no longer requires the C++ extension at import time

Status on this workspace after the transition PR:

- satisfied: 1, 2, 3, 4, and 5 for the current stable NVIDIA runtime
- not yet satisfied: 6; the public package and the `_pyopencl` implementation still import C++ wrapper-owned types
- still unresolved for broader support: the Intel CPU runtime remains unstable for the `localmax` rebuild path

## Post-Switch Follow-Up Worth Keeping In Scope

- improve runtime-selection diagnostics so platform and device enumeration are easier to reconcile across `clinfo`, the legacy wrapper, and PyOpenCL
- separate supported-runtime policy from one-machine-specific observations once broader validation data exists
- revisit whether the Intel CPU runtime instability is a backend bug, a driver issue, or both before treating it as supported