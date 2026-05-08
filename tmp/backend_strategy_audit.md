# Backend Strategy Audit

## Bottom line

clODE no longer has a meaningfully generic backend architecture. It has a single internal OpenCL execution layer implemented through PyOpenCL. Preserving a faux-generic backend seam no longer buys the repo anything.

## Current reality

- `clode.runtime` no longer exposes backend selection at all.
- `clode._opencl` is the canonical internal execution package.
- `clode._opencl/runtime.py` is explicitly single-device and owns a single context and command queue.
- `clode._opencl/source_builder.py`, `registry.py`, and the kernel tree all assume OpenCL C assets and PyOpenCL host objects.
- The public API still carries some broader-looking selection surface, such as `device_ids`, but the implementation does not actually provide general multi-device execution.

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

## Why one might want another backend

- Hardware coverage gaps, especially if OpenCL availability or quality becomes unacceptable on a target platform.
- A CPU reference path for debugging or deterministic validation.
- Vendor-specific performance work that PyOpenCL plus OpenCL cannot reasonably cover.

## Realism and recommendation

The realistic choices are:

1. Stay PyOpenCL-only and simplify the repo around that decision.
2. Add a scoped reference executor later if correctness or testing pressure justifies it.

The unrealistic choice is keeping a generic backend architecture in place today because a future CUDA, SYCL, or Metal backend might someday happen. The current seam is too thin to make that future cheap, and too abstract-looking to clarify the present code.

## Recommended repo decision

- Treat clODE as a PyOpenCL-backed OpenCL package, not a backend-agnostic package.
- Keep `_opencl/` as the internal execution package rather than restoring any generic backend seam.
- Do not preserve interface complexity just to keep hypothetical backend optionality alive.
- Either retire `device_ids` or document it explicitly as compatibility surface rather than real multi-device execution.
- If a reference executor becomes desirable, introduce it as an explicitly scoped correctness tool.

## Consequence for module layout

If the project decides there is no near-term second backend plan, it can reorganize the package around semantic layers instead of pretending the backend seam is the main architecture. That should make the internals easier to read and easier to evolve.
