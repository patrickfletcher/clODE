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
        "trajectory.py",
    ],
    visibility = ["//visibility:public"],
    deps = ["//clode/cpp:py_clode"],
)

py_library(
    name = "clode_all_sources",
    visibility = ["//visibility:public"],
    deps = [
        #"//clode:clode_init",
        "//clode/cpp:py_clode",
        "//clode/python:clode_init",
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
#
#py_wheel(
#    name = "pyclode",
#    # Package data. We're building "example_minimal_package-0.0.1-py3-none-any.whl"
#    distribution = "clode",
#    python_tag = "py3",
#    stamp = 0,
#    strip_path_prefixes = ["src"],
#    version = "0.0.3",
#    deps = [":clode_pkg"],
#)
