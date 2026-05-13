# Backend Strategy Audit

Purpose: preserve the historical rationale for retiring the generic backend seam and standardizing on `_opencl` plus PyOpenCL.
Read when: you need design archaeology for why clODE no longer presents itself as backend-agnostic.
Update when: a signpost to the current source of truth needs correction.

Historical only. The current source of truth is `.design/package_state.md`.

## Bottom line

clODE no longer has a meaningfully generic backend architecture. It has a single internal OpenCL execution layer implemented through PyOpenCL, and the migration-era protocol/factory facade has now been removed.

## Current reality at the time of archival

- `clode.runtime` no longer exposes backend selection at all.
- `clode._opencl` is the canonical internal execution package.
- The public simulators now instantiate the concrete `_opencl` executors directly instead of routing through generic backend protocols and a backend factory shim.
- `clode._opencl/runtime.py` is explicitly single-device and owns a single context and command queue.
- `clode._opencl/source_builder.py`, `registry.py`, and the kernel tree all assume OpenCL C assets and PyOpenCL host objects.
- The public runtime selection surface is explicitly single-device. Any future multi-device work would need a new API and execution model.

## What another backend would actually mean

### 1. Another OpenCL host backend

Example: a direct C API wrapper, cffi layer, Rust extension, or a revived custom extension.

What it would reuse:

- Most of the OpenCL kernels.
- Much of the build-key, source assembly, and problem-shape logic.

What it would still need:

- A full runtime layer for device discovery, contexts, queues, buffers, program builds, and transfers.
- A second struct-layout and kernel-launch path to keep correct and tested.

Assessment:

- Low user value.
- High maintenance duplication.
- Hard to justify when PyOpenCL already exists and is the current dependency story.

### 2. A non-OpenCL accelerator backend

Example: CUDA, HIP, SYCL, Metal, or a JAX/XLA execution path.

What it would reuse:

- Public solver semantics at most.
- Some high-level ideas around problem shapes, observers, and build keys.

What it would still need:

- New kernel sources or a new code-generation target.
- New runtime, memory, and build systems.
- New portability and correctness work for struct layout, math, RNG behavior, and diagnostics.
- Likely a separate optimization and benchmarking effort.

Assessment:

- Potentially valuable if a concrete platform gap appears, especially Apple or CUDA-only demand.
- Realistically a new project, not a small extension of the current seam.

### 3. A CPU or reference backend

Example: NumPy, Numba, or a deliberately slow reference executor.

What it would reuse:

- More of the public API and test references than the two options above.

What it would still need:

- A clear execution-state model.
- A deliberately scoped semantics contract.
- Care around keeping it a correctness tool instead of a second production engine.

Assessment:

- This is the most realistic second implementation if one is ever needed.
- Its value is debugging, correctness, and testing, not high-performance portability.
- It should be framed as a reference executor, not as a justification for preserving a faux-generic production backend layer.

## Historical recommendation

The realistic choices were:

1. Stay PyOpenCL-only and simplify the repo around that decision.
2. Add a scoped reference executor later if correctness or testing pressure justifies it.

The unrealistic choice was keeping a generic backend architecture in place for hypothetical future CUDA, SYCL, or Metal work.
