# PyOpenCL Struct Handling Audit

## Scope

This audit checks whether the current PyOpenCL groundwork is handling OpenCL struct-backed arguments correctly against the current kernel tree and PyOpenCL guidance.

Reference used:

- https://documen.tician.de/pyopencl/howto.html#how-to-use-struct-types-with-pyopencl

## Summary

The transient-side PyOpenCL groundwork was close, but not fully correct in the way PyOpenCL recommends.

The important distinction is:

- `SolverParams` and `ObserverParams` are host-populated structs, so the Python side should use device-matched NumPy dtypes rather than relying on handwritten `align=True` layouts alone.
- `ObserverData` is an observer-specific device-resident struct whose exact layout depends on compile-time macros such as `N_VAR`, `N_AUX`, `N_STORE_EVENTS`, precision, and the selected observer implementation. It should not be treated as a single fixed dtype.

## Findings

### 1. `SolverParams` packing needed a PyOpenCL-native layout match

Status:

- Fixed in the PyOpenCL layer

Before this audit, `clode/_pyopencl/buffers.py` packed `SolverParams` using a handwritten NumPy structured dtype with `align=True`.

That matched the current test device for both single and double precision, but it was still relying on NumPy's host-side ABI assumptions rather than PyOpenCL's device-specific layout matching step.

PyOpenCL's recommended flow is:

1. define a base NumPy dtype
2. call `pyopencl.tools.match_dtype_to_c_struct(...)`
3. use the matched dtype for host packing

The PyOpenCL layer now does exactly that through `clode/_pyopencl/structs.py`.

### 2. `ObserverParams` should follow the same matched-dtype rule

Status:

- Fixed in the PyOpenCL layer as preparation for the feature backend

`ObserverParams` is not yet used by the PyOpenCL execution path, but it is already a host-populated struct passed into kernels as `__constant struct ObserverParams *opars`.

Its field order is stable in the kernel tree:

- four `uint` fields first
- then the precision-dependent `realtype` fields

This one also matches the current test device under `align=True`, but it should still be matched through PyOpenCL rather than assumed.

The PyOpenCL layer now has `get_observer_params_struct(...)` and `pack_observer_params(...)` ready for that later backend work.

### 3. `ObserverData` is the real subtle point, and the old host pattern must not be copied

Status:

- Not implemented on the PyOpenCL path yet
- Existing C++ host-side sizing pattern is not safe to reuse

`ObserverData` is defined differently by each observer implementation under `clode/cpp/observers/*.clh`, and many of the current host-side size calculations are handwritten formulas such as:

```cpp
oi.observerDataSizeDouble = 7*sizeof(cl_double) + sizeof(cl_uint);
```

That pattern is unsafe because struct tail padding matters in double precision.

Measured on the current test device for the simplest case, `ObserverData_basic`:

- handwritten formula: 60 bytes
- PyOpenCL-matched struct size: 64 bytes

So the C++-style byte-count pattern can under-allocate double-precision observer state.

This is the most important audit result.

For the PyOpenCL feature backend, we should not reuse those size formulas. The correct options are:

1. generate a per-observer matched dtype from an explicit Python layout model
2. or generate an exact opaque byte count from the same authoritative layout model that defines the observer struct

What we should not do is keep using `n_real*sizeof(realtype) + n_int*sizeof(uint)` as if it were the struct size.

### 4. We should not inject duplicate struct declarations into the existing kernel sources

Status:

- Confirmed

PyOpenCL's how-to shows prepending the matched C declaration to kernel source. That is correct when the kernel source does not already define the struct.

In clODE, the current source tree already defines:

- `struct SolverParams` in `clode/cpp/clODE_struct_defs.cl`
- `struct ObserverParams` in `clode/cpp/observers.cl`
- observer-specific `ObserverData` typedefs in `clode/cpp/observers/*.clh`

So for this codebase, the right use of `match_dtype_to_c_struct(...)` is host-side dtype matching only. We should keep the existing kernel declarations as the source of truth and avoid prepending duplicate typedefs with the same names.

## Current Rule Set

For the PyOpenCL backend work going forward:

- use `match_dtype_to_c_struct(...)` for any host-populated OpenCL struct
- treat `SolverParams` and `ObserverParams` as matched structured dtypes
- do not treat `ObserverData` as one global dtype
- do not reuse the current C++ observer byte-count formulas on the PyOpenCL path
- keep the current OpenCL source files as the struct declaration authority

## Immediate Impact On PR Sequence

This does not block the transient backend work.

Why:

- PR 9 only needs the transient path, which uses `SolverParams` but not `ObserverParams` or `ObserverData`

It does materially affect the future feature backend work.

Guardrail for later PRs:

- PR 11 should build an explicit observer-struct layout model rather than inheriting `observerDataSize` formulas from the old C++ host layer.