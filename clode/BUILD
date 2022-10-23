load("@rules_python//python:packaging.bzl", "py_package", "py_wheel")

py_library(
    name = "clode_init",
    srcs = [
        "__init__.py",
    ],
    visibility = ["//visibility:public"],
)

py_wheel(
    name = "pyclode",
    # Package data. We're building "example_minimal_package-0.0.1-py3-none-any.whl"
    distribution = "clode",
    python_tag = "py3",
    stamp = 0,
    #strip_path_prefixes = ["src"],
    version = "0.0.8",
    deps = [
        "//clode/python:clode_pkg",
    ],
)
