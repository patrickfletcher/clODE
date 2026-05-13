# PyOpenCL Leverage Audit

Purpose: audit which PyOpenCL services clODE is already using and which remaining host-side helpers are still worth adopting.
Read when: touching runtime selection, program builds, buffer ownership, OpenCL diagnostics, or PyOpenCL-related backlog items.
Update when: clODE adopts a new PyOpenCL helper family or a PyOpenCL leverage backlog item materially changes.

## Bottom line

The live code is already using PyOpenCL at the right level for the solver hot path: explicit `Context`, `CommandQueue`, `Program`, `Kernel`, `Buffer`, `enqueue_copy`, and `enqueue_nd_range_kernel` control remains appropriate for clODE's custom kernels.

The audit found only a narrow helper layer beyond those raw primitives:

- `pyopencl.tools.match_dtype_to_c_struct(...)` in `clode/_opencl/structs.py`, `test/test_opencl_structs.py`, and `tools/probe_opencl_runtime.py`
- `pyopencl.characterize.has_src_build_cache(...)` in `clode/_opencl/runtime.py`
- raw build-log extraction through `Program.get_build_info(..., program_build_info.LOG)` in `clode/_opencl/program_cache.py`

No repo code currently uses `cache_dir`, `MemoryPool`, `ImmediateAllocator`, `capture_call`, `enqueue_fill_buffer`, `enqueue_fill`, `enqueue_map_buffer`, `pyopencl.array`, elementwise/reduction/scan helpers, `clrandom`, `pyopencl.cltypes`, or PyOpenCL's pytest parametrization helpers.

## Where clODE is already using PyOpenCL well

### 1. Raw runtime and kernel control

- `clode/_opencl/runtime.py` creates `Context` and `CommandQueue` explicitly.
- `clode/_opencl/program_cache.py` builds `Program` objects and materializes named `Kernel` handles.
- `clode/_opencl/buffers.py` owns raw `Buffer` allocation and host/device copies.
- `clode/_opencl/executors.py` launches kernels directly with `enqueue_nd_range_kernel(...)`.

This remains the right design boundary. clODE's solver, observer, and source-assembly semantics are package-specific and should stay above PyOpenCL rather than being forced into array-side helper abstractions.

### 2. Struct-layout matching

- `clode/_opencl/structs.py` uses `match_dtype_to_c_struct(...)` to keep Python-side struct packing aligned with the selected device.
- `test/test_opencl_structs.py` verifies the matched dtype sizes and packed values.
- `tools/probe_opencl_runtime.py` uses the same helper for diagnostic probes.

This is already the correct PyOpenCL leverage point for the current struct path. There is no current repo evidence that `pyopencl.cltypes` would replace it cleanly or reduce maintenance.

### 3. Basic diagnostics and inventory

- `clode/runtime/query.py` uses raw `get_platforms()` / `get_devices()` and device attributes for public inventory reporting.
- `clode/_opencl/program_cache.py` surfaces compiler/build logs through `get_build_info(...)`.
- `tools/probe_opencl_runtime.py` provides a PyOpenCL-first diagnostic path for runtime and struct-layout failures.

### 4. Limited `characterize` usage

- `clode/_opencl/runtime.py` already logs whether the selected device/runtime reports source-build cache support via `has_src_build_cache(...)`.
- No broader `pyopencl.characterize` helpers are used yet for tuning or diagnostics.

## High-confidence leverage gaps

### 1. Build-cache controls and diagnostics

`clode/_opencl/program_cache.py` currently calls `program.build(options=...)` without threading `cache_dir`. clODE's `BuildKey` and `SourceBundle` logic should remain, but PyOpenCL's own source-build cache controls and diagnostics should be the first stop before growing more custom build-cache behavior.

### 2. Buffer clearing and transfer helpers

`clode/_opencl/buffers.py` still clears observer data by allocating host zeros and uploading them with `enqueue_copy(...)`. No repo code uses `enqueue_fill_buffer`, `enqueue_fill`, or `enqueue_map_buffer`. This is the clearest concrete place where PyOpenCL host-side helpers might reduce unnecessary host traffic.

### 3. Memory pooling

`clode/_opencl/buffers.py` still allocates every buffer manually. No repo code uses `ImmediateAllocator`, `MemoryPool`, `SVMAllocator`, or `SVMPool`. Pooling is worth benchmarking for transient, trajectory, and feature buffer churn, but only if profiling shows real allocation pressure.

### 4. Kernel repro capture

No repo code uses `Kernel.capture_call()`. This looks most useful as a developer-only debug path for vendor/runtime failures rather than as part of the default runtime surface.

## Lower-confidence or lower-priority leverage areas

- `pyopencl.array`, `ElementwiseKernel`, reductions, scans, and `clrandom` are not used anywhere in the repo.
- The current execution path is still better served by raw custom kernels than by array-side PyOpenCL helpers.
- If these helpers are adopted, the strongest current fit is diagnostics, preprocessing, tests, or narrow support utilities, not the solver hot path.
- PyOpenCL's pytest parametrization helpers are also unused today; they are relevant only if clODE intentionally broadens its internal device-matrix testing.

## Recommendation

Keep raw `Program` / `Kernel` / `Buffer` / queue control as the runtime model. Narrow active PyOpenCL leverage work to four concrete gaps:

1. build-cache controls and related diagnostics (`cache_dir`, broader `characterize` helpers)
2. device-side fill/map helpers instead of host-zero uploads where it matters
3. optional memory-pool experiments for measured buffer churn
4. optional `capture_call()` support for debugging hard runtime failures

Treat array-side helpers, external RNG helpers, and device-matrix pytest helpers as secondary audits for diagnostics/testing/support workflows unless the live code develops a stronger need.
