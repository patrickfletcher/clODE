---
applyTo: "clode/kernels/**"
description: Use when editing OpenCL kernels under clode/kernels/. Read the single-precision ODE solver guide first and keep kernel numerics lean, stable, and evidence-backed.
---

# Kernel Numerics Rules

- Read `.design/reference/single_precision_ode_solver_guide.md` before changing solver math, reductions, tolerances, scaling, time integration, or mixed-precision behavior in `clode/kernels/**`.
- For repo-specific float32 behavior, time bookkeeping, observer numerics, interpolation, or feature accumulation, also read `.design/reference/single_precision_numerics_note.md` when relevant.
- Let each stepper or observer own only the state and stage-time geometry it actually needs; keep shared wrappers minimal.
- Avoid numerically hostile algebra, naive reductions, and tolerance tightening that is not justified in float32.
- When kernel numerics behavior or guidance changes, update the smallest relevant `.design/reference/` note in the same PR.