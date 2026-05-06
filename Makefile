PYTHON ?= python
DIST_DIR ?= dist

.PHONY: install-test install-docs install-dev test-smoke test-opencl test-extended test-long build-docs build-dist legacy-cpp-wrapper paper

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

test-opencl:
	$(PYTHON) tools/run_test_bundle.py opencl

test-extended:
	$(PYTHON) tools/run_test_bundle.py extended

test-long:
	$(PYTHON) tools/run_test_bundle.py long

build-docs:
	$(PYTHON) -m mkdocs build --strict

build-dist:
	$(PYTHON) -m build --outdir $(DIST_DIR)
	$(PYTHON) -m twine check $(DIST_DIR)/*

legacy-cpp-wrapper:
	$(MAKE) -C clode/cpp install-wrapper PYTHON=$(PYTHON)

paper:
	docker run --rm \
    --volume ./paper:/data \
    --user $(id -u):$(id -g) \
    --env JOURNAL=joss \
    openjournals/inara
