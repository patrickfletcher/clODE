cc_library(
    name = "cl",
    hdrs = ["cl.hpp"],
)

cc_library(
    name = "cl2",
    hdrs = ["cl2.hpp"],
)

cc_library(
    name = "openclhpp",
    hdrs = ["opencl.hpp"],
)

config_setting(
    name = "macos_opencl",
    constraint_values = [
        "@platforms//os:macos",
    ],
    visibility = ["//visibility:public"],
)

cc_library(
    name = "OpenCL",
    includes = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["external/com_github_khronosgroup_opencl"],
    }),
    linkopts = select({
        ":macos_opencl": ["-framework OpenCL"],
        "//conditions:default": [],
    }),
    visibility = ["//visibility:public"],
    deps = [
        ":cl",
        ":cl2",
        ":openclhpp",
    ] + select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["@com_github_khronosgroup_opencl_windows//:opencl_windows"],
    }),
)
