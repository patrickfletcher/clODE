load("@rules_python//python:packaging.bzl", "py_package", "py_wheel")

py_wheel(
    name = "pyclode",
    # Package data. We're building "example_minimal_package-0.0.1-py3-none-any.whl"
    distribution = "clode",
    python_tag = "py3",
    stamp = 0,
    version = "0.2.0",
    deps = [
        "//clode/python:clode_pkg",
    ],
)
