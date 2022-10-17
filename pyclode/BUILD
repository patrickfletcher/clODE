load("@rules_python//python:packaging.bzl", "py_package", "py_wheel")

py_library(
    name = "clode_init",
    srcs = [
        "__init__.py",
        "features.py",
        "observer.py",
        "problem_info.py",
        "runtime.py",
        "solver_params.py",
        "stepper.py",
    ],
    visibility = ["//visibility:public"],
    deps = ["//src:py_clode"],
)

py_package(
    name = "clode_pkg",
    # Only include these Python packages.
    #packages = ["examples.wheel"],
    deps = [
        ":clode_init",
        "//src:py_clode",
    ],
)

py_wheel(
    name = "pyclode",
    # Package data. We're building "example_minimal_package-0.0.1-py3-none-any.whl"
    distribution = "clode",
    python_tag = "py3",
    version = "0.0.2",
    deps = [":clode_pkg"],
)
