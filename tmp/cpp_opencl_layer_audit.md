# C++ OpenCL Layer Audit

## Scope

This report covers the current C++ and OpenCL solver layer under `clode/cpp` and the immediate Python wrappers in `clode/runtime.py`, `clode/solver.py`, `clode/trajectory.py`, and `clode/features.py`.

The goal is to document how the current system is assembled, identify what is working well, identify what is driving maintenance cost, and recommend what should and should not survive a move to PyOpenCL.

## Executive Summary

The current design is conceptually sound at the numerical level. The OpenCL assets are organized around a small number of kernel entry points, a reusable stepper library, a reusable observer library, and a clean ensemble execution model where one work-item advances one ODE instance. That core model is worth preserving.

The maintenance cost is concentrated in the host layer, not in the numerical kernels themselves. The current C++ layer owns device selection, program compilation, buffer allocation, host/device marshaling, and Python bindings, but most of that logic is mechanical and duplicated again in the Python wrappers. It also relies heavily on compile-time defines, mixed-language metadata files, hand-maintained size formulas, and rebuild rules that are only partially enforced by the Python API.

The recommended migration is not a kernel rewrite first. It is a host-runtime rewrite first. Keep the OpenCL kernel assets and the high-level Python API shape, move runtime ownership into Python with PyOpenCL, replace the C++ manifest and buffer-management layer with a Python program cache plus buffer manager, and only then simplify the kernel source layout where that clearly reduces coupling.

## Current Structure

### Main host-side components

| Component | Role | Key files |
| --- | --- | --- |
| `OpenCLResource` | OpenCL platform/device discovery, context creation, queue creation, program build, device info queries | `clode/cpp/OpenCLResource.hpp`, `clode/cpp/OpenCLResource.cpp` |
| `CLODE` | Base simulator for transient solves, base buffer ownership, build options, kernel construction, host/device state transfer | `clode/cpp/CLODE.hpp`, `clode/cpp/CLODE.cpp` |
| `CLODEtrajectory` | Extends `CLODE` with trajectory storage buffers and the `trajectory` kernel | `clode/cpp/CLODEtrajectory.hpp`, `clode/cpp/CLODEtrajectory.cpp` |
| `CLODEfeatures` | Extends `CLODE` with observer state, feature buffers, observer initialization, and the `features` kernel | `clode/cpp/CLODEfeatures.hpp`, `clode/cpp/CLODEfeatures.cpp` |
| `clode_cpp_wrapper` | pybind11 bridge exporting C++ structs, runtime selection, and simulator classes to Python | `clode/cpp/CLODEpython.cpp` |
| Python wrappers | Public user API, ensemble shaping, array reshaping, RHS conversion, runtime initialization | `clode/runtime.py`, `clode/solver.py`, `clode/trajectory.py`, `clode/features.py` |

### Main kernel-side components

| Component | Role | Key files |
| --- | --- | --- |
| Entry-point kernels | Top-level kernels for transient, trajectory, observer warmup, and features | `clode/cpp/transient.cl`, `clode/cpp/trajectory.cl`, `clode/cpp/initializeObserver.cl`, `clode/cpp/features.cl` |
| Common types and helpers | Precision typedefs, solver parameter structs, interpolation and running-stat helpers, RNG | `clode/cpp/realtype.cl`, `clode/cpp/clODE_struct_defs.cl`, `clode/cpp/clODE_utilities.cl`, `clode/cpp/clODE_random.cl` |
| Stepper library | Method registry plus per-method implementation headers and shared wrappers | `clode/cpp/steppers.cl`, `clode/cpp/steppers/*.clh` |
| Observer library | Observer registry plus per-observer state and logic headers | `clode/cpp/observers.cl`, `clode/cpp/observers/*.clh` |
| User RHS | Problem-specific `getRHS(...)` function loaded from a `.cl` file or generated from Python/XPP | runtime-selected source file, `clode/function_converter.py`, `clode/xpp_parser.py` |

## How The Current System Fits Together

### 1. Problem definition and runtime selection

The Python API constructs a `ProblemInfo` object containing:

- the path to the OpenCL RHS source file
- `nVar`, `nPar`, `nAux`, and `nWiener`
- variable, parameter, and auxiliary names

The runtime is selected through `OpenCLResource`, either by platform/device id or by vendor/device type filters.

### 2. Base program assembly in C++

`CLODE` loads `transient.cl` into `clprogramstring` during construction. `CLODEtrajectory` appends `trajectory.cl`. `CLODEfeatures` appends `initializeObserver.cl` and `features.cl`.

Separately, `setProblemInfo` reads the user RHS source file into `ODEsystemsource`.

At build time, the final program source is:

```text
clprogramstring + ODEsystemsource
```

That means the top-level kernel entry-point files are concatenated directly, while the common kernel library is pulled in through OpenCL `#include` directives.

### 3. Build options drive specialization

`CLODE::setCLbuildOpts` generates compile flags for:

- precision: `CLODE_SINGLE_PRECISION` or `CLODE_DOUBLE_PRECISION`
- stepper selection: for example `EXPLICIT_RK4` or `EXPLICIT_DOPRI5`
- problem dimensions: `N_VAR`, `N_PAR`, `N_AUX`, `N_WIENER`
- include path: `-I<clodeRoot>`

`CLODEfeatures::buildCL` adds:

- observer selection define, for example `USE_OBSERVER_LOCAL_MAX`
- `N_STORE_EVENTS` for event timestamp storage

This compile-time specialization is central to the design. It is how the kernels get fixed-size private arrays such as `realtype xi[N_VAR]`, fixed-size observer structs, and selected stepper and observer implementations.

### 4. Include-based kernel composition

The composition chain for the core kernels looks like this:

```text
transient.cl
  -> clODE_random.cl
  -> clODE_struct_defs.cl
  -> clODE_utilities.cl
  -> realtype.cl
  -> steppers.cl
       -> selected stepper .clh
       -> shared fixed/adaptive stepper wrapper
  -> user getRHS(...) from the problem source

trajectory.cl
  -> same common includes as transient.cl

initializeObserver.cl
  -> same common includes
  -> observers.cl
       -> all observer headers are visible
       -> only the selected observer defines the active ObserverData typedef and functions

features.cl
  -> same common includes
  -> observers.cl
```

Two details matter here.

First, `steppers.cl` and `observers.cl` have a dual role. Under `__cplusplus` they expose host-side metadata helpers. Under OpenCL compilation they act as include hubs for the kernel implementations. That is clever, but it tightly couples host metadata and kernel source layout.

Second, user RHS code is not included through the same library structure. It is appended as a separate trailing source string and is expected to define `getRHS(...)` with the exact signature the steppers call.

## Host And Device State Model

### Ensemble execution model

The execution model is straightforward and good.

- one global work-item corresponds to one independent ODE solve
- ensemble size is `nPts`
- each work-item copies its parameters, initial state, and RNG state into private arrays
- the solver loop then runs entirely in private memory except for final write-back and optional trajectory or feature output

This is the core abstraction worth preserving.

### Data layout conventions

The host Python wrappers flatten arrays in Fortran order. The kernels then treat data as variable-major arrays:

- state and parameters: `buffer[var_index * nPts + point_index]`
- trajectory samples: `buffer[store_index * nPts * width + item_offset * nPts + point_index]`
- feature output: `F[feature_index * nPts + point_index]`

This layout is internally consistent and favors coalesced access across ensemble members. It should remain an internal implementation detail, not something users need to understand.

### Base buffers

`CLODE` owns the common buffers:

- `d_tspan`
- `d_sp`
- `d_pars`
- `d_x0`
- `d_xf`
- `d_RNGstate`
- `d_dt`
- `d_tf`

The host mirrors are stored as `std::vector<double>` and cast down to `float` for device transfer when single precision is enabled.

### Trajectory specialization

`CLODEtrajectory` adds:

- `d_t`
- `d_x`
- `d_dx`
- `d_aux`
- `d_nStored`

The trajectory kernel stores the initial point and then every `nout`th accepted step up to `max_store`.

### Feature specialization

`CLODEfeatures` adds:

- `d_odata` for persistent observer state
- `d_op` for observer parameters
- `d_F` for feature output

The observer lifecycle is:

1. allocate `ObserverData` storage per ensemble member
2. optionally run `initializeObserver` for two-pass detectors
3. run `features`
4. update `ObserverData` in global memory so later calls can continue without reinitializing unless requested

That continuation model is valuable and should be kept.

## What The Current Design Does Well

- The numerical decomposition is strong. Steppers, observer logic, and kernel entry points are separated in a way that makes the solver logic readable once the build model is understood.
- The ensemble execution model is simple and efficient. One work-item per solve is the right mental model for this problem class.
- Compile-time specialization gives the kernels fixed-size private arrays and observer structs, which keeps the hot loops simple and efficient.
- The observer contract is well-defined. Each observer provides the same conceptual hooks: initialize, optional warmup, per-step update, event detection, event feature extraction, finalization.
- The Python-side RHS generation tools are high-value assets. `function_converter.py` and `xpp_parser.py` already move problem definition into Python, which is exactly the direction a PyOpenCL port should continue.
- The code exposes the final built program string, which is useful when debugging generated source and build options.
- The public Python API is already higher level than the C++ layer. That means the C++ layer can be removed without forcing a user-facing API reset.

## What Is Driving Maintenance Cost

### 1. The host layer is doing mostly mechanical work

The C++ layer mainly owns OpenCL context creation, buffer allocation, data copies, build option assembly, and pybind exposure. None of that is the scientific core. It is the part most likely to become easier to maintain in Python with PyOpenCL.

### 2. Host metadata and kernel source are coupled together

`steppers.cl` and `observers.cl` act as both host metadata registries and OpenCL include hubs. That creates a hidden dependency between:

- C++ string maps of available methods
- feature-name generation
- observer-data size formulas
- actual OpenCL struct definitions and function bodies

This saves duplication in one place and creates it in another. The observer size calculations, in particular, are hand-maintained formulas that must stay aligned with the actual `ObserverData` struct layout.

### 3. Compile-time specialization is useful, but the rebuild surface is broad

A rebuild is required when any of these change:

- precision
- stepper
- observer type
- problem dimensions
- `maxEventTimestamps`

That is manageable, but only if the host runtime owns it explicitly. Right now the rules are spread across C++ comments, Python flags, and wrapper behavior.

### 4. Python invalidation is only partially wired through

The Python wrappers track `_cl_program_is_valid`, but the lazy rebuild path in `transient`, `trajectory`, and `features` is commented out. That means the API contains the concept of rebuild invalidation, but the enforcement is incomplete.

For a migration, this is a signal that rebuild logic should become a first-class concern of the runtime rather than a side flag.

### 5. Several capabilities are broader on paper than in practice

- `OpenCLResource` can create a context with multiple devices, but the solver classes launch kernels through `getQueue()` with the default device index and do not partition work across devices.
- Trajectory storage is preallocated and not chunked, even though the source already contains TODOs acknowledging that long trajectories should be streamed or paged.
- Observer timestamp storage is embedded into `ObserverData`, which makes storage capacity part of the compile contract rather than a runtime output policy.

### 6. Boilerplate is duplicated across kernels

`transient.cl`, `trajectory.cl`, `initializeObserver.cl`, and `features.cl` each repeat the same initial loading pattern:

- copy `tspan`, parameters, state, RNG state into private memory
- compute the initial RHS
- run a time loop
- write back final state, RNG state, `dt`, and final time

The duplication is not catastrophic, but it raises the cost of changing the shared execution model.

### 7. Capacity and storage semantics are harder to reason about than they should be

Two examples stand out.

- Trajectory storage uses `max_store` as both a solver/storage control and an allocation size, while the kernel also stores the initial point. That makes the true capacity contract easy to misread and difficult to evolve safely.
- Observer event storage size is a compile-time property because event arrays live inside `ObserverData`. That couples analysis detail to kernel recompilation.

### 8. The build and packaging stack is heavier than the runtime problem requires

The current Python package build uses Bazel plus pybind11 to expose a runtime that is, at its core, an OpenCL program builder, buffer manager, and kernel launcher. That is a large amount of tooling around a workload PyOpenCL already supports directly.

## What Is Worth Porting

The following parts are worth preserving almost directly.

- The ensemble execution model.
- The public Python simulator classes and their user-facing semantics.
- The kernel math library: `realtype.cl`, `clODE_random.cl`, `clODE_utilities.cl`.
- The stepper implementations in `clode/cpp/steppers/*.clh`.
- The observer implementations in `clode/cpp/observers/*.clh`, at least as the phase-one behavioral reference.
- The continuation semantics for `x0`, `dt`, final time, RNG state, and observer state.
- The Python AST-to-OpenCL RHS conversion path.
- The ability to inspect or emit the final generated source for debugging.

## What Should Be Simplified In A PyOpenCL Port

### 1. Replace the C++ runtime with a Python runtime

The new runtime should own:

- context and queue creation
- device queries
- source assembly
- build-option generation
- program caching
- buffer allocation and reuse
- host/device transfers
- build log reporting

That removes `OpenCLResource`, `CLODE`, `CLODEtrajectory`, `CLODEfeatures`, and the pybind module from the critical path.

### 2. Move stepper and observer metadata into Python

Instead of keeping host metadata behind `__cplusplus` branches in `.cl` files, create Python registries such as:

- available steppers and the define each one requires
- available observers
- observer feature-name builders
- observer build keys
- observer state sizing rules or explicit layouts

The OpenCL source tree should then contain only OpenCL concerns.

### 3. Keep compile-time specialization in phase one

A full move to runtime-generic dimensions is not the best first target. The current kernels depend heavily on compile-time sizes for private arrays and observer state. Replacing that with runtime-generic kernels would change both performance characteristics and kernel structure.

The practical migration path is:

- keep compile-time specialization for precision, stepper, observer, and problem dimensions
- make that specialization explicit in a Python program-cache key
- cache programs by something like `(precision, stepper, observer, n_var, n_par, n_aux, n_wiener, n_store_events)`

That preserves the current numerical structure while removing the C++ scaffolding around it.

### 4. Separate solver state from optional output capacity

Two simplifications are especially valuable.

- For trajectories, separate integration control from output storage capacity. Use chunked host readback or paged device buffers instead of forcing one monolithic `max_store` allocation.
- For observers, move optional event timestamp storage out of `ObserverData` when practical. Small persistent observer state and larger optional event-output buffers should not be the same concern.

This change would reduce rebuild pressure and make long-run analysis much easier to support.

### 5. Centralize rebuild decisions

The runtime should have one explicit rule: if the build key changes, rebuild or fetch from cache; otherwise reuse the compiled program.

That is cleaner than the current distributed logic where rebuild-sensitive fields exist in C++, some invalidation flags exist in Python, and lazy rebuild hooks are partially disabled.

### 6. Keep the memory layout internal and well documented

The current column-major flattening is reasonable for the current kernel access pattern. It does not need to be user-visible. Preserve the layout internally if benchmarking continues to support it, but make the reshape and flatten rules live in one Python buffer module rather than across wrappers and C++ methods.

### 7. Start with explicit single-device support

The current code does not truly implement multi-device execution, even though the runtime can create multi-device contexts. The PyOpenCL port should either:

- explicitly support one device at a time in phase one, or
- implement real partitioned execution across devices

The current halfway state is not worth carrying forward.

## What Should Probably Be Dropped

- The pybind11 wrapper layer.
- The Bazel-based extension build as a required runtime dependency for the Python API.
- Mixed host/kernel manifest files that depend on `__cplusplus` branching.
- The nominal multi-device API unless real device partitioning is implemented.
- Duplicated default values for solver and observer parameters across language boundaries.

## Recommended PyOpenCL Architecture

### Phase-one target

Keep the user-facing Python API roughly as it is today and replace the backend with a small set of Python runtime objects.

Suggested structure:

| Python component | Responsibility |
| --- | --- |
| `runtime.py` | Context and queue creation, device selection, build-log handling |
| `kernel_registry.py` | Stepper and observer metadata, source fragments, build-key rules |
| `program_cache.py` | Build-key to `pyopencl.Program` cache |
| `buffers.py` | Allocation, resizing, flattening, reshaping, host/device transfer helpers |
| `simulator.py` | Base simulator orchestration |
| `trajectory.py` | Trajectory-specific output management |
| `features.py` | Observer lifecycle and feature output management |

### Source assembly strategy

The most streamlined source strategy is:

1. keep the existing OpenCL kernel library files largely intact for phase one
2. remove host-only metadata from those files
3. let Python assemble the top-level program source and build options
4. cache the resulting program by build key

That preserves the battle-tested numerical code while making the assembly path explicit and inspectable in Python.

### Migration order

1. Reproduce the current build key and kernel assembly in Python with PyOpenCL.
2. Reuse the current kernel assets unchanged as much as possible.
3. Match current simulator outputs against the C++ backend.
4. Once parity exists, simplify the source tree and output/storage model.
5. Only then consider deeper kernel redesigns such as runtime-generic dimensions or merged entry-point kernels.

## Bottom Line

The current C++ layer is not the part to preserve. The current OpenCL solver model is.

Port the kernel library and execution model, not the C++ scaffolding. Keep compile-time specialization where it still buys clarity and performance, but move ownership of build keys, source assembly, buffer management, and rebuild policy into Python. Simplify metadata, separate persistent solver state from optional output storage, and treat true multi-device support as either an explicit feature or out of scope.

That path removes the largest maintenance burden while preserving the numerical behavior and most of the code that is actually unique to clODE.