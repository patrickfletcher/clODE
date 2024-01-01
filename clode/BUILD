load("@rules_python//python:packaging.bzl", "py_package", "py_wheel")

py_library(
    name = "clode_init",
    srcs = [
        "__init__.py",
        "features.py",
        "function_converter.py",
        "runtime.py",
        "solver.py",
        "trajectory.py",
        "xpp_parser.py",
    ],
    visibility = ["//visibility:public"],
    deps = ["//clode/cpp:py_clode"],
)

py_library(
    name = "clode_all_sources",
    visibility = ["//visibility:public"],
    deps = [
        ":clode_init",
        "//clode/cpp:py_clode",
    ],
)

py_package(
    name = "clode_pkg",
    visibility = ["//visibility:public"],
    # Only include these Python packages.
    #packages = ["examples.wheel"],
    deps = [
        ":clode_all_sources",
    ],
)

py_wheel(
    name = "pyclode",
    # Package data. We're building "example_minimal_package-0.0.1-py3-none-any.whl"
    distribution = "clode",
    python_tag = "py3",
    stamp = 0,
    version = "",
    deps = [
        ":clode_pkg",
    ],
)
