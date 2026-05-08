# Legacy Surface Cleanup Plan

Status: active forward-looking cleanup plan for the first repository change set that removes the legacy C++ wrapper from mainline clODE.

## Scope and assumptions

- PyOpenCL rollout is complete and is the default product path.
- This plan assumes `_CLODE_BACKEND=cpp` will be removed from the main repository rather than promoted to a long-term supported extra.
- Generated build output under `bazel-*` and `site/` is local housekeeping, not part of the source cleanup itself.

## Tmp cleanup completed on 2026-05-07

Archived today:

- `tmp/archived/backend_migration_history_2026_05_07/backend_core_test_suite.md`
- `tmp/archived/backend_migration_history_2026_05_07/backend_overhaul_scope_map.md`

Archive location:

- `tmp/archived/backend_migration_history_2026_05_07/`

Keep in `tmp/` for now:

- `tmp/backend_pr_task_breakdown.md`
- `tmp/pyopencl_rollout_guardrails.md`
- `tmp/pyopencl_packaging_release_audit.md`
- `tmp/ci_cleanup_plan.md`
- `tmp/pyopencl_backend_design.md`
- `tmp/pyopencl_post_migration_plan.md`
- `tmp/cpp_opencl_layer_audit.md`
- `tmp/pyopencl_struct_audit.md`
- `tmp/python_simulation_flow_reference.md`

The last three are still useful while the comparison path exists. Once the wrapper-removal PR lands, they should move to `tmp/archived/` as well.

## Wrapper-coupled surfaces that must change before removal

| Area | Paths | Why it blocks removal | Planned action |
| --- | --- | --- | --- |
| Backend selection and wrapper loading | `clode/runtime.py`, `clode/_backends/factory.py` | Runtime still accepts `cpp`, checks for wrapper binaries, and lazy-loads the wrapper module | Remove the `cpp` backend name and wrapper-loader helpers; simplify factory code to a single PyOpenCL path |
| C++ backend adapter | `clode/_backends/cpp.py` | Entire module is a wrapper adapter around `clode.cpp.clode_cpp_wrapper` | Delete the module |
| Wrapper type bridges | `clode/types.py` | `*_to_cpp` and `*_from_cpp` helpers exist only to translate through the wrapper boundary | Delete the bridge helpers and keep the Python-owned dataclasses |
| Comparison build/install tooling | `tools/install_cpp_wrapper.py`, `clode/cpp/Makefile` | These commands exist only to build and drop the wrapper into a source checkout | Delete them if the wrapper leaves the main repo |
| Comparison tests and markers | `test/test_pyopencl_transient_backend.py`, `test/test_pyopencl_trajectory_backend.py`, `test/test_pyopencl_feature_backend.py`, `test/conftest.py`, `tools/run_test_bundle.py`, `pyproject.toml`, `Makefile`, `test/README.md` | The repo still treats wrapper comparison as a supported validation mode | Remove the three comparison tests, the `legacy_cpp_comparison` marker, the bundle, and the related contributor docs |
| Published compatibility docs | `docs/install.md`, `docs/init_runtime.md`, `docs/querying_opencl.md` | User docs still describe how to build and compare against the wrapper | Remove wrapper build instructions and wrapper-specific runtime-ordering notes |

## Pure legacy or decision-gated surfaces

| Bucket | Paths | Current role | Recommended action once the wrapper is dropped |
| --- | --- | --- | --- |
| Bazel workspace and dependency glue | `WORKSPACE`, `bazel/**`, `bazelisk.py`, `clode/cpp/BUILD`, `matlab/BUILD`, `samples/BUILD` | Retained only for wrapper, standalone C++, and MEX builds | Delete if the repository no longer ships the old native stack |
| C++ host runtime and wrapper sources | `clode/cpp/**` | Legacy OpenCL host runtime, pybind11 wrapper, and standalone C++ build surface | Preferred: remove from the main repo entirely; fallback: split to a separate compatibility repository/package |
| MATLAB/MEX interface | `matlab/**`, `docs/matlab.md`, `mkdocs.yml`, `.github/ISSUE_TEMPLATE/bug_report.md` | Unmaintained interface built on the old C++ path | Delete or quarantine with the wrapper; remove docs nav and stale support prompts |
| C++ and MEX sample tree | `samples/**` | Legacy C++ and MATLAB sample harnesses plus old output references | Delete after checking whether any `.cl` or `.ode` assets should be rehomed into `examples/` or `docs/` |
| Generated local outputs | `bazel-bin/`, `bazel-out/`, `bazel-clODE/`, `bazel-testlogs/`, `site/` | Local build artifacts; not part of the maintained source tree | Keep ignored and delete locally as needed; do not carry them as migration work |
| Paper pipeline | `paper/**`, `.github/workflows/draft_pdf.yml` | Research-paper artifact pipeline, not part of the wrapper path | Keep unless the paper workflow itself is being retired |

## Recommended cleanup sequence

### Phase 1: remove the comparison backend from Python

- Remove `cpp` backend selection from `clode/runtime.py` and `clode/_backends/factory.py`.
- Delete `clode/_backends/cpp.py`.
- Delete the C++ bridge helpers from `clode/types.py`.
- Update runtime error text and public docs so PyOpenCL is the only in-tree backend story.

### Phase 2: remove comparison-only tests and tooling

- Delete `test/test_pyopencl_transient_backend.py`.
- Delete `test/test_pyopencl_trajectory_backend.py`.
- Delete `test/test_pyopencl_feature_backend.py`.
- Remove the `legacy_cpp_comparison` marker, bundle, Makefile target, and related test documentation.
- Delete `tools/install_cpp_wrapper.py` and wrapper-install references from docs.

### Phase 3: decide the fate of the old native tree

Preferred path:

- Delete `clode/cpp/**`, `WORKSPACE`, `bazel/**`, and `bazelisk.py` from the main repository.

Fallback path if there is real downstream demand:

- Move the old native stack to a separate compatibility repository or package instead of leaving it half-maintained inside the default repo.

### Phase 4: remove orthogonal legacy interfaces that depend on the old native tree

- Delete `matlab/**`.
- Delete `docs/matlab.md` and remove the mkdocs nav entry.
- Remove the stale Matlab prompt from `.github/ISSUE_TEMPLATE/bug_report.md`.
- Delete `samples/**` unless specific reusable assets are intentionally relocated first.

### Phase 5: finish repository simplification

- Remove stale compatibility language from docs and contributor notes.
- Collapse test documentation around a single backend.
- Archive the remaining wrapper-era tmp notes after the deletion PR lands.
- Leave `site/` and Bazel output cleanup as local housekeeping rather than source edits.

## Concrete delete list once the drop is approved

Delete immediately with wrapper removal:

- `clode/_backends/cpp.py`
- the `cpp` branch in `clode/_backends/factory.py`
- wrapper-loader helpers and `cpp` backend acceptance in `clode/runtime.py`
- `problem_info_to_cpp`, `solver_params_to_cpp`, `observer_params_to_cpp`, and the matching `*_from_cpp` helpers in `clode/types.py`
- `tools/install_cpp_wrapper.py`
- `clode/cpp/Makefile`
- `test/test_pyopencl_transient_backend.py`
- `test/test_pyopencl_trajectory_backend.py`
- `test/test_pyopencl_feature_backend.py`
- the `legacy_cpp_comparison` marker and bundle references in `pyproject.toml`, `test/conftest.py`, `tools/run_test_bundle.py`, `Makefile`, and `test/README.md`
- wrapper-comparison sections in `docs/install.md`, `docs/init_runtime.md`, and `docs/querying_opencl.md`

Delete if the main repository will no longer host any old native compatibility surface:

- `WORKSPACE`
- `bazel/**`
- `bazelisk.py`
- `clode/cpp/**`
- `matlab/**`
- `samples/**`
- `docs/matlab.md`
- the Matlab nav entry in `mkdocs.yml`
- the Matlab prompt in `.github/ISSUE_TEMPLATE/bug_report.md`

Keep unless a separate retirement decision is made:

- `clode/kernels/**`
- `examples/**`
- `paper/**`
- `.github/workflows/draft_pdf.yml`

## Pre-delete checks

- Verify that no import of `clode.cpp` or `clode_cpp_wrapper` remains in `clode/`, `test/`, `tools/`, or `docs/`.
- Decide whether any files under `samples/` should move into `examples/` or `docs/` before deletion.
- Decide whether Matlab support is being removed outright or archived elsewhere with the wrapper.
- Decide whether the old native stack needs a separate compatibility repository before deleting Bazel and `clode/cpp/`.

## Success criteria for the cleanup PR

- `resolve_backend_name()` no longer accepts `cpp`.
- No `clode.cpp` or `clode_cpp_wrapper` import remains in the maintained package path.
- No `legacy_cpp_comparison` marker or bundle remains.
- Published docs describe one backend only.
- If the native compatibility stack is dropped, no Bazel workspace files remain in the main repo outside archived history notes.
