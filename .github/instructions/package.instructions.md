---
applyTo: "clode/**"
description: Use when editing clODE package source or kernels. Preserve the canonical package layout and compatibility strategy.
---

# Package Layout Rules

- Prefer canonical homes: `clode.problem`, `clode.observers`, `clode.simulation`, `clode.runtime`, and `clode._opencl`.
- Treat flat root modules such as `clode.solver`, `clode.features`, `clode.trajectory`, `clode.types`, `clode.function_converter`, `clode.xpp_parser`, and `clode.opencl_builtins` as compatibility barrels unless the task is explicitly about them.
- `_opencl` owns runtime, build, buffer, cache, and dispatch mechanics; semantic definitions should live on the Python side when a clear home exists.
- `clode/kernels/` is shipped runtime code. Keep kernel changes aligned with build keys, observer definitions, and stepper definitions.
- When editing `clode/kernels/**`, read `.design/reference/single_precision_ode_solver_guide.md` first and then the smallest repo-specific numerics note needed for the task.
- Keep the public API stable unless the task clearly requires a change, and update the relevant `.design` doc when semantic boundaries move.