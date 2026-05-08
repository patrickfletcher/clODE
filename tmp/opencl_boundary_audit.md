# OpenCL Boundary Audit

## Scope

Audit the post-refactor boundary between the canonical internal execution layer in `_opencl/` and the public semantic subpackages `problem/`, `simulation/`, `observers/`, and `runtime/`.

This audit is intentionally narrower than the earlier module-layout plan. The question here is not whether `_opencl/` should exist; it should. The question is whether any remaining contents in `_opencl/` actually belong in a public semantic package, or whether any public packages are still carrying implementation-era concerns.

## Completed State

- Public semantic ownership now lives in `clode.problem`, `clode.simulation`, `clode.observers`, and `clode.runtime`.
- Internal execution ownership now lives in `clode._opencl`.
- `clode._pyopencl` is now only a compatibility package-level shim.
- `_backends/` is no longer used as a live simulation/runtime seam; only compatibility wrappers remain.

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

## Real Remaining Candidates

### 1. `clode._opencl.observer_metadata`

Status:

- this is the strongest remaining boundary question

Why it is mixed:

- part of the file is clearly observer-domain logic: per-observer feature naming, count logic, and observer-specific data-shape decisions
- part of the file is clearly execution-layer logic: matching observer data structs against OpenCL layouts for a specific runtime and precision

Recommendation:

- do not move the whole file into `observers/`
- instead, split it in a future pass:
  - pure observer catalog/name/schema logic into `clode.observers` (likely `definitions.py` or a new `metadata.py`)
  - runtime/precision-specific struct matching stays in `_opencl/`

Reason:

- moving the whole module outward would leak OpenCL struct-layout concerns into the public observer namespace
- keeping the whole module inward hides observer semantics that are not actually execution-specific

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

- package ownership is now semantic, but several internal symbol names still reflect the previous package name

Examples:

- `PyOpenCLTransientBackend`
- `PyOpenCLFeatureBackend`
- `PYOPENCL_BACKEND_VERSION`
- `PyOpenCLDependencyError`

Recommendation:

- no urgent correctness issue
- migrate callers and tests toward the new canonical `_opencl` aliases first, then decide whether the old symbol names should be retired or kept indefinitely as compatibility aliases

Reason:

- this is naming debt, not architectural debt
- removing it immediately would add churn without improving the module boundary itself

### 4. Public re-export of `_load_pyopencl` and `_require_pyopencl`

Status:

- these remain reachable from `clode.runtime` via `runtime/__init__.py`

Recommendation:

- likely safe to keep for now because they are underscore-prefixed and already clearly internal
- later, consider stopping the re-export from `clode.runtime.__init__` if there is no real downstream use

Reason:

- the functions are binding-specific implementation helpers, not part of the semantic public runtime surface

### 5. Compatibility packages still present

Status:

- `clode._pyopencl` remains as a package-level shim
- `_backends/` remains for compatibility, primarily `rhs.py` and historical factory/protocol imports

Recommendation:

- keep for now
- remove only after internal tests/docs and any downstream imports stop relying on them

Reason:

- these are migration tails, not semantic placement mistakes

## Non-Candidates

These items may look close to public semantics, but should stay internal:

- `_opencl.structs`: depends on runtime device layout matching and precision-specific C struct packing
- `_opencl.source_builder`: kernel assembly and compile defines are internal execution concerns
- `_opencl.program_cache`: purely internal build/runtime caching
- `_opencl.buffers`: purely internal device memory ownership and transfer layout

## Suggested Follow-up Order

1. If another semantic pass is desired, split `observer_metadata.py` rather than moving it wholesale.
2. After compatibility appetite is clearer, decide whether `_pyopencl` and `_backends` shims can be removed.
3. Only after that, decide whether historical `PyOpenCL...` symbol names should be retired in favor of `_opencl` aliases.

## Validation Snapshot

Validation run after this refactor:

- targeted pytest bundle covering backend-rhs, runtime, `_opencl` internals, and simulation/backend contracts: `41 passed`
- smoke bundle: `21 passed`