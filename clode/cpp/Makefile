PYTHON ?= python
BAZEL ?= bazel
REPO_ROOT := $(abspath $(CURDIR)/../..)

.PHONY: build-wrapper install-wrapper

build-wrapper:
	cd $(REPO_ROOT) && $(BAZEL) build //clode/cpp:clode_cpp_wrapper

install-wrapper: build-wrapper
	cd $(REPO_ROOT) && $(PYTHON) tools/install_cpp_wrapper.py
