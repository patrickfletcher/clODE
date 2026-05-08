# OpenCL Boundary Audit

## Scope

Audit the post-refactor boundary between the canonical internal execution layer in `_opencl/` and the public semantic subpackages `problem/`, `simulation/`, `observers/`, and `runtime/`.

This audit is intentionally narrower than the earlier module-layout plan. The question here is not whether `_opencl/` should exist; it should. The question is whether any remaining contents in `_opencl/` actually belong in a public semantic package, or whether any public packages are still carrying implementation-era concerns.

## Completed State

- Public semantic ownership now lives in `clode.problem`, `clode.simulation`, `clode.observers`, and `clode.runtime`.
- Internal execution ownership now lives in `clode._opencl`.
- The observer feature-name and two-pass catalog logic now lives in `clode.observers.metadata`.
- The old binding-named compatibility package and `_backends/` have been removed.
- The runtime package no longer re-exports or depends on backend-selection compatibility helpers.
- Historical `PyOpenCL...` symbol names inside `_opencl` have been retired in favor of `OpenCL...` names.

## Confirmed Good Boundaries

These modules look correctly placed and should stay internal to `_opencl/`:

- `buffers.py`: device buffers and host-device copy layout are execution-layer concerns.
- `runtime.py`: PyOpenCL context/queue/device selection and cache ownership are execution-layer concerns.
- `program_cache.py`: build/program caching is internal compile/runtime machinery.
- `source_builder.py`: build-option synthesis and kernel-source assembly are internal compile machinery.
- `structs.py`: OpenCL/C struct matching and packed transfer layouts are execution details, even though they consume public dataclasses.
- `models.py`: `BuildKey`, `SourceBundle`, `ProgramBundle`, and `ProblemShape` are internal compile/cache models rather than public problem or simulation concepts.

These public modules also look correctly scoped after the refactor:

- `problem/*`: problem definition, source loading, XPP/Python ingestion, and authoring helpers.
- `simulation/*`: simulator orchestration, solver params, and result containers.
- `observers/*`: observer enum and observer parameter schema.
- `runtime/*`: public runtime/device selection, querying, and logging.

## Remaining Candidates

### 1. `clode._opencl.observer_metadata`

Status:

- the strongest boundary question from the previous audit has now been partially resolved

Why it is mixed:

- part of the file is clearly observer-domain logic: per-observer feature naming, count logic, and observer-specific data-shape decisions
- part of the file is clearly execution-layer logic: matching observer data structs against OpenCL layouts for a specific runtime and precision

Completed change:

- pure observer catalog/name/two-pass logic has been moved into `clode.observers.metadata`
- runtime/precision-specific struct matching remains in `_opencl/observer_metadata.py`

Reason:

- this keeps OpenCL struct-layout concerns out of the public observer namespace without hiding the observer catalog itself

### 2. `clode._opencl.registry`

Status:

- acceptable where it is, but semantically mixed

Why it is mixed:

- it contains compile-time stepper and observer macro mappings, which reference simulation and observer concepts
- but those mappings are specifically about OpenCL kernel entrypoints and compile defines

Recommendation:

- keep it in `_opencl/` for now
- revisit only if the observer catalog or stepper catalog needs a public or shared declarative source of truth

Reason:

- the current data is tightly coupled to kernel filenames and compile defines, so it still reads as internal execution metadata rather than public API metadata

### 3. Historical `PyOpenCL...` names inside `_opencl`

Status:

- resolved in the current follow-up pass

### 4. Public re-export of binding-loader helpers

Status:

- resolved in the current follow-up pass; these helpers are no longer re-exported from `clode.runtime`

### 5. Compatibility packages still present

Status:

- resolved in the current follow-up pass; both compatibility packages have been removed

## Non-Candidates

These items may look close to public semantics, but should stay internal:

- `_opencl.structs`: depends on runtime device layout matching and precision-specific C struct packing
- `_opencl.source_builder`: kernel assembly and compile defines are internal execution concerns
- `_opencl.program_cache`: purely internal build/runtime caching
- `_opencl.buffers`: purely internal device memory ownership and transfer layout

## Suggested Follow-up Order

1. Revisit `clode._opencl.registry` only if stepper or observer catalog data needs a shared declarative source of truth.
2. If observer design work continues, separate public observer definitions further from runtime struct-layout concerns.

## Validation Snapshot

Validation run after this refactor:

- targeted pytest bundle covering backend-rhs, runtime, `_opencl` internals, and simulation/backend contracts: `41 passed`
- smoke bundle: `22 passed`