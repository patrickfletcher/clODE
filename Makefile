PYTHON ?= python
DIST_DIR ?= dist

.PHONY: install-test install-docs install-dev test-smoke test-frontend test-runtime-api test-numerics test-release test-opencl test-opencl-internal test-extended build-docs build-dist paper

install-test:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e ".[test]"

install-docs:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e ".[docs]"

install-dev:
	$(PYTHON) -m pip install --upgrade pip
	$(PYTHON) -m pip install -e ".[dev]"

test-smoke:
	$(PYTHON) tools/run_test_bundle.py smoke

test-frontend:
	$(PYTHON) tools/run_test_bundle.py frontend

test-runtime-api:
	$(PYTHON) tools/run_test_bundle.py runtime_api

test-numerics:
	$(PYTHON) tools/run_test_bundle.py numerics

test-release:
	$(PYTHON) tools/run_test_bundle.py release

test-opencl:
	$(PYTHON) tools/run_test_bundle.py opencl

test-opencl-internal:
	$(PYTHON) tools/run_test_bundle.py opencl_internal

test-extended:
	$(PYTHON) tools/run_test_bundle.py extended

build-docs:
	$(PYTHON) -m mkdocs build --strict

build-dist:
	$(PYTHON) -m build --outdir $(DIST_DIR)
	$(PYTHON) -m twine check $(DIST_DIR)/*

paper:
	docker run --rm \
    --volume ./paper:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/inara
