# PyOpenCL Leverage Audit

## Bottom line

clODE is using PyOpenCL at the right level for its custom kernels, but there are several host-side services that PyOpenCL already provides and clODE should prefer before growing more custom runtime machinery.

## Where clODE is already using PyOpenCL well

- Raw `Program`, `Kernel`, `Buffer`, `Context`, and `CommandQueue` usage is appropriate for the custom solver, trajectory, and observer kernels.
- `pyopencl.tools.match_dtype_to_c_struct` is already being used correctly for device-matched struct layout in `clode/_opencl/structs.py`.
- Direct kernel control keeps the performance-critical integration path in clODE's hands, which PyOpenCL is meant to enable rather than replace.

## High-value leverage points

### 1. Build caching

PyOpenCL already caches built binaries on disk when `Program.build()` is called on source programs, and also supports an explicit `cache_dir`. clODE's `ProgramCache` currently provides only an in-memory `BuildKey -> ProgramBundle` cache.

Recommendation:

- Keep the `BuildKey` and `SourceBundle` logic, because that is clODE-specific.
- Do not grow `ProgramCache` into a second full compiler-cache system.
- Make PyOpenCL's disk cache behavior part of clODE's runtime diagnostics and developer workflow.

### 2. Memory pools

`clode/_opencl/buffers.py` currently performs raw buffer allocation and manual reuse decisions. PyOpenCL already provides `ImmediateAllocator`, `MemoryPool`, `SVMAllocator`, and `SVMPool`.

Recommendation:

- Prototype `MemoryPool(ImmediateAllocator(queue))` for common transient, feature, and trajectory buffer churn.
- Keep clODE's shape-aware buffer ownership model, but stop reinventing allocation policy if pooling shows clear wins.

### 3. Device characterization

PyOpenCL exposes `pyopencl.characterize` helpers for build-cache detection, SIMD group estimates, local-memory limits, and fast-math option discovery.

Recommendation:

- Use these helpers in diagnostics and tuning experiments before building more custom device heuristics.
- This is especially useful for future work on local-memory usage, work-group tuning, and fast-math tradeoffs.

### 4. Debug reproduction helpers

`Kernel.capture_call()` can emit a self-contained PyOpenCL repro for a failing kernel invocation.

Recommendation:

- Add a developer-facing path for capturing minimal repros when a vendor runtime or observer path misbehaves.
- This is a good fit for hard-to-reproduce OpenCL build or execution bugs.

### 5. Testing and device parametrization

PyOpenCL ships a pytest parametrization helper that can help fan tests out across devices and platforms.

Recommendation:

- Consider this for signal-only internal test lanes.
- Keep the current narrow release gate selection, but use PyOpenCL's test helpers when intentionally widening device coverage.

### 6. Transfer and clearing helpers

clODE currently clears some device buffers by allocating host zeros and uploading them. PyOpenCL already provides `enqueue_fill_buffer`, `enqueue_fill`, `enqueue_map_buffer`, and richer transfer helpers.

Recommendation:

- Prefer device-side fills for large zero-initialization paths such as observer-state clears.
- Use mapping helpers when debugging or when a host-visible buffer workflow is clearer than a copy-heavy path.

### 7. RNG and array-side utilities

PyOpenCL exposes array utilities, elementwise kernels, reductions, scans, and `clrandom` generators backed by Random123 components.

Recommendation:

- Treat these as support tools for tests, preprocessing, diagnostics, and smaller helper operations.
- Do not assume they are drop-in replacements for clODE's per-step solver-integrated RNG path.

## Places where clODE still risks reinventing the wheel

- `clode/_opencl/program_cache.py` should stay a small bundle cache, not become a parallel build-cache subsystem.
- `clode/_opencl/buffers.py` should not drift into a homemade memory-pool framework if PyOpenCL's pool APIs are sufficient.
- Device diagnostics should use `pyopencl.characterize` where possible instead of growing ad hoc capability logic.
- Large host-zero uploads should be rechecked against `enqueue_fill_buffer`.

## Where custom clODE code is still justified

- Source assembly and build keys remain clODE-specific because they encode solver, observer, precision, and problem-shape semantics.
- Kernel registries and executors remain clODE-specific because PyOpenCL is not supposed to know clODE's numerical model.
- Struct packing helpers still belong in clODE even though the low-level layout matching comes from PyOpenCL.

## Recommendation

Use PyOpenCL more aggressively for host-side services and diagnostics, but keep clODE's custom kernels and execution semantics firmly in clODE. The main design rule should be: if the logic is OpenCL-host plumbing, first ask whether PyOpenCL already provides it.
