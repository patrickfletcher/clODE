# Contributing

Contributions are welcome. Keep changes scoped, keep public wording aligned with the current Python package, and update tests or docs when user-visible behavior changes.

## Development setup

clODE currently targets Python 3.10 and newer.

For a full development environment from a source checkout:

```bash
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
python -m pip check
```

The repository also provides equivalent Make targets:

```bash
make install-dev
make install-docs
make install-test
```

An OpenCL runtime is required for runtime-backed tests and examples. Before debugging clODE itself, verify what the host exposes:

```bash
python -c "import clode; print(clode.query_opencl())"
python tools/probe_opencl_runtime.py
```

If you do not have a working OpenCL runtime on the current machine, stick to the smoke bundle and docs build.

## Common commands

```bash
make test-smoke
make test-frontend
make test-runtime-api
make test-numerics
make test-release
make test-opencl-internal
make build-docs
```

The authoritative test bundle definitions live in `tools/run_test_bundle.py`. Useful direct invocations are:

```bash
python tools/run_test_bundle.py smoke
python tools/run_test_bundle.py release
python tools/run_test_bundle.py opencl_internal
python -m mkdocs build --strict
```

Bundle intent:

- `smoke`: driver-independent packaging and frontend checks.
- `frontend`: user-facing conversion and ingestion checks.
- `runtime_api`: public runtime and simulator contracts.
- `numerics`: numerical regression coverage.
- `release`: the OpenCL-backed release gate.
- `opencl_internal`: focused tests for internal OpenCL support layers.

## Docs and public surfaces

Keep the public-facing surfaces distinct:

- `README.md`: GitHub landing page and package index long description.
- `docs/index.md`: documentation home page.
- `paper/paper.md`: research-facing narrative.

When a change affects user workflows, examples, or project positioning, update the relevant public surface in the same pass. Keep public wording centered on the current Python package and supported workflows; do not add historical implementation archaeology to user-facing pages. Because `README.md` also serves as the package index long description, keep it concise and prefer absolute URLs over repo-relative links.

## Maintainer planning docs

The repo-local planning surface lives under `.design/`.

- Start with `.design/README.md` and stop at the smallest read set that fits the task.
- Keep the root `.design/` docs small and authoritative; move deep dives to `.design/reference/`, history to `.design/archived/`, and scratch work to `.design/tmp/`.
- If a change affects package layout, the active implementation target, or the meaning of a reference note, update the smallest authoritative `.design` doc in the same pass.

## Pull requests

- Include tests when behavior changes or regressions are fixed.
- Update documentation when public APIs, workflows, or installation guidance change.
- Prefer canonical package modules for new implementation work rather than compatibility barrels.
- Call out the platform, device, and runtime used for OpenCL-backed changes when it matters for reproduction.

For bug reports and reproduction notes, include at least:

- operating system
- Python version
- clODE version or commit
- OpenCL platform, device, and driver/runtime details
- minimal script or model needed to reproduce the issue

## Scope notes

Citation metadata and the paper are still in progress. Do not add provisional citation text or release-facing scholarly claims unless the supporting paper and metadata are ready at the same time.
