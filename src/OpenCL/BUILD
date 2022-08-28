cc_library(
    name = "cl",
    hdrs = ["cl.hpp"],
)

cc_library(
    name = "cl2",
    hdrs = ["cl2.hpp"],
)

config_setting(
    name = "macos_opencl",
    constraint_values = [
        "@bazel_tools//platforms:osx",
    ],
    visibility = ["//visibility:public"],
)

config_setting(
    name = "linux_opencl",
    constraint_values = [
        "@bazel_tools//platforms:linux",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "OpenCL",
    linkopts = select({
        ":linux_opencl": [],
        ":macos_opencl": ["-framework OpenCL"],
        "//conditions:default": [],
    }),
    visibility = ["//visibility:public"],
    deps = [
        ":cl",
        ":cl2",
    ],
    #    deps = select({
    #        ":linux_opencl": [
    #            ":cl",
    #            "cl2",
    #        ],
    #        ":macos_opencl": [],
    #        "//conditions:default": [
    #            ":cl",
    #            "cl2",
    #        ],
    #    }),
)
