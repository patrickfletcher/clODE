# PyOpenCL Post-Migration Plan

## Purpose

This note is about the state after the migration is no longer merely "PyOpenCL works" and is instead "the project is structurally organized around a Python-owned PyOpenCL backend".

It answers four questions:

1. what benefits we actually want from relying on PyOpenCL
2. how close the current codebase is to those benefits
3. what the post-migration codebase and packaging should look like
4. what to simplify or defer once the legacy C++ wrapper is no longer the default path

## The Benefits We Actually Want

The point of adopting PyOpenCL is not just to replace one runtime API with another. The real goals are:

### 1. A normal Python install story

- no required Bazel build for the default package
- no required pybind11 wrapper in the default install path
- a pure-Python wheel plus packaged OpenCL assets
- a release process that behaves like a normal modern Python project

### 2. Python ownership of runtime behavior

- explicit device selection, source assembly, build keys, and caches in Python
- easier debugging of program text, build options, and runtime selection
- less logic hidden behind a compiled wrapper

### 3. Better modeling of solver state

- clearer separation between problem definition, solver configuration, runtime state, and cached outputs
- fewer duplicated defaults across language boundaries
- less ad hoc invalidation logic spread across multiple layers

### 4. Better modeling of observers

- observer metadata and state treated as first-class Python-owned concepts
- a clear distinction between persistent observer state and optional event-output storage
- a viable path to user-customizable observers without forcing full trajectory storage

### 5. Lower maintenance cost

- fewer build systems and language boundaries
- simpler CI and release workflows
- a smaller set of runtime-critical files to reason about during debugging

## Where We Are Now

The current state is strong on backend parity and noticeably weaker on packaging and post-migration cleanup.

### Already landed

- the internal backend seam exists
- the PyOpenCL backend runs transient, trajectory, and feature paths
- the public Python types and runtime facade are now Python-owned
- `import clode` no longer requires the C++ extension when the PyOpenCL path is selected
- the 76-test PyOpenCL acceptance bundle is green on the stable local NVIDIA runtime

### Not landed yet

- the default package is still Bazel-backed
- the default wheel is still not pure Python
- kernel assets still live under `clode/cpp/`
- the default backend has not switched to PyOpenCL
- user-authored observers are not a polished public API
- solver state and observer state are still only partially modeled as explicit Python concepts

## How Close Are We?

The answer depends on which finish line matters.

### PyOpenCL backend parity

Close.

The core numerical and behavioral migration work is substantially done. The current system already proves that a Python-owned PyOpenCL backend can execute the existing solver and observer kernel paths with parity-grade tests.

### PyOpenCL as the default shipped product

Not there yet, but the remaining work is narrower and better defined.

The blockers are now mostly packaging, asset layout, release workflow, and runtime-support policy. That is a much better problem to have than unresolved numerical parity.

### Post-migration streamlined architecture

Still meaningfully incomplete.

The codebase still carries transition structure, legacy packaging assumptions, and some API surfaces that exist mainly because the wrapper used to own them.

### Easy custom observers and richer solver-state modeling

Not close yet as a user-facing feature.

The ingredients are present:

- Python-owned observer metadata
- Python-owned runtime and build models
- Python-owned config structs

But the public abstraction for authoring or loading custom observers does not exist yet.

## Recommended Sequence From Here

### Step 1: Finish PR15 cleanly

This is still the highest-leverage next move.

PR15 should:

- remove Bazel from the default package build path
- move kernel assets into explicit package data
- stop shipping the compiled extension and C++ source tree in the default wheel
- clean up version sourcing and release workflows

This is the step that turns the project from "PyOpenCL-capable" into "PyOpenCL-owned".

### Step 2: Switch the default backend only after packaging is clean

Only then should the default backend move to PyOpenCL.

Switching the backend first would improve the default runtime choice but would still leave the package structurally tied to the legacy extension. That is the wrong order.

### Step 3: Do a focused post-migration cleanup pass

After the default switch lands and stabilizes, then simplify:

- remove or isolate the legacy C++ backend
- simplify the internal package layout
- reduce transition-only compatibility logic
- tighten docs around one supported default path

### Step 4: Start the modeling and customization phase

Only after the package and runtime path are simplified should the project take on bigger design work around solver state and custom observers.

That is where the real long-term value lies, but it is easier to do after the runtime ownership story is fully settled.

## Recommended Post-Migration Folder Structure

There are two credible options.

## Option A: Minimal-diff cleanup

Keep most current names and only remove legacy pieces.

```text
clode/
  __init__.py
  types.py
  runtime.py
  solver.py
  trajectory.py
  features.py
  function_converter.py
  xpp_parser.py
  opencl_builtins.py
  _backends/
    protocol.py
    factory.py
    rhs.py
  _pyopencl/
    runtime.py
    registry.py
    source_builder.py
    program_cache.py
    buffers.py
    structs.py
    observer_metadata.py
    executors.py
    models.py
    errors.py
  kernels/
    *.cl
    steppers/*.clh
    observers/*.clh
```

Pros:

- minimal churn after a risky migration
- low rename cost
- easier diff against the transition codebase

Cons:

- `_pyopencl` stays as a transition-era name even after PyOpenCL is the only backend
- internal build/runtime/source concerns remain spread across a package whose name describes the implementation, not the role

## Option B: Cleaner long-term internal layout

Rename the implementation packages around responsibilities rather than backend history.

```text
clode/
  __init__.py
  types.py
  runtime.py
  solver.py
  trajectory.py
  features.py
  function_converter.py
  xpp_parser.py
  opencl_builtins.py
  kernels/
    *.cl
    steppers/*.clh
    observers/*.clh
  _internal/
    rhs.py
    registry.py
    source_builder.py
    runtime.py
    program_cache.py
    buffers.py
    structs.py
    observer_metadata.py
    executors.py
    models.py
    errors.py
```

Pros:

- more honest final architecture
- clearer ownership boundaries for runtime, build, and observer machinery
- easier to reason about once the C++ path is gone

Cons:

- rename churn across many imports
- little immediate value before the package and default backend are stable

## Recommendation

Take Option A through PR15 and the default-backend switch, then decide whether Option B still feels worth the rename churn.

That is the pragmatic choice. The codebase does not need a large rename while packaging and rollout are still in flight.

## What Packaging Should Look Like

The packaging target should be a conventional Python package.

## Recommended end state

- `pyproject.toml` is the authoritative packaging configuration
- the default wheel is `py3-none-any`
- runtime kernel assets are packaged explicitly as package data
- `pyopencl` is a normal runtime dependency for the default package once PyOpenCL is the default backend
- the OpenCL driver and ICD remain an environment prerequisite, documented clearly but not bundled
- release publishing happens only from a dedicated tag-driven workflow
- versioning comes from one authoritative source

## Recommendation on `pyopencl` dependency shape

Once PyOpenCL is the default backend, `pyopencl` should be a required dependency of the main package, not an optional extra.

Pros:

- the default install path is honest about what is required to run the package
- fewer branches in docs and CI
- avoids the confusing state where the default package installs but the default backend dependency is absent

Cons:

- users on systems without a usable OpenCL runtime will still need environment setup
- some environments may still find PyOpenCL installation or runtime setup non-trivial

Recommendation:

- make `pyopencl` required when the default backend flips
- keep the runtime-environment troubleshooting in the docs, not in package extras

## Recommendation on legacy C++ compatibility

If a legacy C++ fallback is still worth keeping for one or two releases, it should not remain welded into the default package build.

Better options are:

- a separate compatibility package
- an explicit legacy extra with clearly isolated build logic
- a separate branch if no active users need packaged legacy support

The cleanest long-term answer is a separate compatibility package. Keeping the legacy backend inside the default build path defeats most of the packaging win.

## Better Modeling of Solver Components

The current public types are a start, not the full model.

## Recommended internal models

### Problem definition

Keep `ProblemInfo`, but consider eventually splitting it into:

- `ProblemDefinition`: names, dimensions, noise count, RHS source metadata
- `EnsembleData`: host-side variable and parameter arrays plus ensemble shape

Why:

- today, problem metadata and live ensemble values are still orchestrated together inside simulator logic
- splitting them makes cache keys, rebuild decisions, and state reuse easier to reason about

### Solver configuration versus solver state

`SolverParams` should remain the configuration object, but the runtime should also have an explicit internal state object, something like:

- current `t_span`
- live `dt` buffer state
- RNG state residency
- build key and compiled-program identity
- cached final-state or trajectory buffers

Why:

- the current continuation semantics are valuable, but some of them still live as implicit device-side behavior plus cache invalidation rules
- a modeled solver state makes those rules easier to document and extend

### Output policies

Trajectory storage policy should eventually become distinct from integration policy.

Today, `max_store` still couples solver configuration to output allocation. Long term, those should be separated so chunked or paged trajectory storage can exist without overloading core solver configuration.

## Better Modeling of Observers

This is the largest post-migration design opportunity.

## What is good about the current observer concept

- low memory footprint
- continuation-friendly device state
- much better scaling than storing full trajectories for many analyses
- a natural fit for OpenCL execution

## What is still missing

- a clean public definition of an observer
- a stable Python-owned manifest for feature names, persistent state layout, and optional event outputs
- a separation between persistent observer state and optional event timestamp storage
- a user-facing extension path for custom observers

## Recommended staged design

### Stage A: internal observer definition model

Create an internal `ObserverDefinition` concept owned in Python, with fields like:

- observer name
- kernel fragments or source entrypoints
- build-time defines
- feature-name builder
- persistent state layout
- optional event-output layout
- whether an initialization pass is required

Pros:

- makes built-in observers explicit and easier to test
- removes more host metadata from kernel source files
- creates the foundation for customization later

Cons:

- requires a careful refactor of existing observer metadata and source assembly

### Stage B: split persistent state from optional event storage

This should happen before exposing custom observers publicly.

Pros:

- reduces rebuild pressure tied to event timestamp retention
- makes observer memory costs easier to reason about
- improves the path to more flexible long-run analyses

Cons:

- non-trivial buffer and kernel interface refactor

### Stage C: experimental custom observer API

Only after Stages A and B are stable should the project expose a low-level experimental custom observer interface.

Recommendation:

- start with a manifest-plus-kernel-fragment API for advanced users
- do not promise a high-level Python DSL too early

That keeps the first customizable observer story honest and debuggable.

## What To Remove If The C++ Wrapper Is Fully Dropped

These are the obvious removal candidates:

- `clode/_backends/cpp.py`
- the pybind extension module and its type stubs
- host-side C++ files under `clode/cpp/` that only exist for the wrapper runtime
- Bazel build files and Bazelisk glue from the default Python build path
- Bazel-centric publish workflows

These should not necessarily disappear immediately on day one:

- the backend protocol and factory
- parity-style tests that compare old and new behavior

Recommendation:

- keep the backend seam until at least one clean post-switch release has shipped
- remove it only when it is clear that a second backend is not coming back and the abstraction no longer pays for itself

## Codebase Streamlining Opportunities After The Switch

Once the default backend and packaging are stable, the next cleanup pass should target these items:

### 1. Move kernel assets out of `clode/cpp/`

This is the most obvious structural mismatch left over from the old architecture.

### 2. Collapse duplicated transition logic

Remove code paths that only exist to preserve the wrapper boundary or support two host runtimes.

### 3. Deprecate or clarify nominal multi-device arguments

Either implement real multi-device execution or document and deprecate the parts of the API that imply it without delivering it.

### 4. Centralize defaults

Make `SolverParams` and `ObserverParams` the single source of truth for defaults and remove any remaining duplicated values.

### 5. Simplify docs around one supported path

The docs should stop telling two different installation and runtime stories once the default backend switch is complete.

## Recommended Near-Term Decisions

1. Do not let post-migration renaming block PR15 packaging cleanup.
2. Move kernel assets into package data as part of PR15, even if the public API does not change.
3. Treat PyOpenCL as the required runtime dependency once the default backend flips.
4. Keep custom observers out of the public API until the internal observer-definition model is explicit.
5. Plan one deliberate post-switch cleanup pass instead of letting transition structure linger indefinitely.

## Bottom Line

The project is already close to the important numerical and runtime milestone: PyOpenCL can act as the real backend.

It is not yet close to the packaging and architectural finish line that would make the migration feel complete.

The next work should therefore stay disciplined:

- finish packaging and asset ownership first
- switch the default backend second
- simplify the codebase third
- only then open the deeper design phase around solver-state modeling and customizable observers

That ordering preserves momentum and gives the project a cleaner base for the more ambitious post-migration improvements.
