cc_library(
    name = "opencl_resources",
    srcs = ["OpenCLResource.cpp"],
    hdrs = ["OpenCLResource.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        "//clode/cpp/OpenCL",
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

cc_library(
    name = "cl_definitions",
    hdrs = [
        "observers.cl",
        "odedriver.cl",
        "realtype.cl",
        "steppers.cl",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "clode_cl_api",
    hdrs = [
        "clODE_random.cl",
        "clODE_struct_defs.cl",
        "clODE_utilities.cl",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":cl_definitions",
        "//clode/cpp/OpenCL",
    ],
)

cc_library(
    name = "clode_api",
    hdrs = [
        "CLODE.hpp",
    ],
    visibility = ["//visibility:public"],
    deps = [
        "opencl_resources",
        ":cl_definitions",
        ":clode_cl_api",
        "//clode/cpp/OpenCL",
    ],
)

cc_library(
    name = "clode",
    srcs = ["CLODE.cpp"],
    copts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [
        ],
        "@platforms//os:windows": ["/DCOMPILING_DLL"],
    }),
    visibility = ["//visibility:public"],
    deps = [
        ":clode_api",
        ":opencl_resources",
        "//clode/cpp/OpenCL",
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

cc_library(
    name = "clode_features",
    srcs = ["CLODEfeatures.cpp"],
    hdrs = ["CLODEfeatures.hpp"],
    copts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [
        ],
        "@platforms//os:windows": ["/DCOMPILING_DLL"],
    }),
    deps = [
        ":clode",
        ":clode_api",
        ":opencl_resources",
        "//clode/cpp/OpenCL",
        "//clode/cpp/observers",
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

cc_library(
    name = "clode_trajectory",
    srcs = ["CLODEtrajectory.cpp"],
    hdrs = ["CLODEtrajectory.hpp"],
    copts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [
        ],
        "@platforms//os:windows": ["/DCOMPILING_DLL"],
    }),
    deps = [
        ":clode",
        ":clode_api",
        ":opencl_resources",
        "//clode/cpp/OpenCL",
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

filegroup(
    name = "opencl",
    srcs = [
        "clODE_random.cl",
        "clODE_struct_defs.cl",
        "clODE_utilities.cl",
        "features.cl",
        "initializeObserver.cl",
        "observers.cl",
        "odedriver.cl",
        "realtype.cl",
        "steppers.cl",
        "trajectory.cl",
        "transient.cl",
        "//clode/cpp/observers:cl",
        "//clode/cpp/steppers:cl",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "cpp",
    includes = ["."],
    visibility = ["//visibility:public"],
    deps = [
        ":clode",
        ":clode_features",
        ":clode_trajectory",
    ],
)

cc_binary(
    name = "clode_cpp_wrapper",
    srcs = ["CLODEpython.cpp"],
    includes = [
        ".",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [
            "-lOpenCL",
        ],
        "@platforms//os:windows": [],
    }),
    linkshared = True,
    linkstatic = select({
        "@platforms//os:osx": False,
        "@platforms//os:linux": True,
        "@platforms//os:windows": True,
    }),
    visibility = ["//visibility:public"],
    deps = [
        ":clode",
        ":clode_features",
        ":clode_trajectory",
        "//clode/cpp/logging:python_sink",
        "@com_github_pybind_pybind11//:pybind11",
        "@local_config_python//:python_headers",
    ],
)

# This wrapper exists to expose the clode.so binary to python
py_library(
    name = "py_clode",
    data = [
        ":clode_cpp_wrapper",
        ":opencl",
        "//samples:lactotroph",
    ],
    imports = ["."],
    visibility = ["//visibility:public"],
)
