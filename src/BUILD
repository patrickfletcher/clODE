cc_library(
    name = "opencl_resources",
    srcs = ["OpenCLResource.cpp"],
    hdrs = ["OpenCLResource.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        "//src/OpenCL",
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
        "//src/OpenCL",
    ],
)

cc_library(
    name = "clode_api",
    hdrs = [
        "CLODE.hpp",
    ],
    visibility = ["//visibility:public"],
    deps = [
        ":cl_definitions",
        ":clode_cl_api",
        "//src/OpenCL",
    ],
)

cc_library(
    name = "clode",
    srcs = ["CLODE.cpp"],
    visibility = ["//visibility:public"],
    deps = [
        ":clode_api",
        ":opencl_resources",
        "//src/OpenCL",
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

cc_library(
    name = "clode_features",
    srcs = ["CLODEfeatures.cpp"],
    hdrs = ["CLODEfeatures.hpp"],
    deps = [
        ":clode",
        ":clode_api",
        ":opencl_resources",
        "//src/OpenCL",
        "//src/observers",
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

cc_library(
    name = "clode_trajectory",
    srcs = ["CLODEtrajectory.cpp"],
    hdrs = ["CLODEtrajectory.hpp"],
    deps = [
        ":clode",
        ":clode_api",
        ":opencl_resources",
        "//src/OpenCL",
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
        "//src/observers:cl",
        "//src/steppers:cl",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "src",
    includes = ["."],
    visibility = ["//visibility:public"],
    deps = [
        ":clode",
        ":clode_features",
        ":clode_trajectory",
    ],
)

#apple_dynamic_framework_import(
#    name = "Python_H",
#    framework_imports = glob(["Python.framework/**"]),
#)

cc_binary(
    name = "clode.so",
    srcs = ["CLODEpython.cpp"],
    #    copts = select({
    #        #"@platforms//os:osx": ["-I/Library/Frameworks/Python.framework/Versions/Current/Headers/"],
    #        "@platforms//os:osx": [
    #            "-F/Library/Frameworks",
    #            "-framework Python",
    #        ],
    #        "@platforms//os:linux": [],
    #        "@platforms//os:windows": [],
    #    }),
    includes = [
        ".",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL", "-static"],
        "@platforms//os:windows": [],
    }),
    linkshared = True,
    visibility = ["//visibility:public"],
    deps = [
        #":Python_H",
        ":clode",
        ":clode_features",
        ":clode_trajectory",
        "//src/logging:python_sink",
        "@com_github_pybind_pybind11//:pybind11",
        "@local_config_python//:python_headers",
        #"@org_python_ftp//:python",
    ],
)

py_library(
    name = "py_clode",
    data = [
        ":clode.so",
        ":opencl",
        "//samples:lactotroph",
    ],
    imports = ["."],
    visibility = ["//visibility:public"],
)
