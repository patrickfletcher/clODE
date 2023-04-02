filegroup(
    name = "lactotroph",
    srcs = ["lactotroph.cl"],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "testTrans",
    size = "large",
    srcs = ["testTransient.cpp"],
    data = [
        "lactotroph.cl",
        "//clode/cpp:cl_definitions",
        "//clode/cpp:clode_cl_api",
        "//clode/cpp:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//clode/cpp",
    ] + select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["@com_github_khronosgroup_opencl//:opencl_headers"],
    }),
)

cc_test(
    name = "testTraj",
    srcs = ["testTrajectory.cpp"],
    data = [
        "lactotroph",
        "//clode/cpp:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//clode/cpp",
    ] + select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["@com_github_khronosgroup_opencl//:opencl_headers"],
    }),
)

cc_test(
    name = "testFeat",
    srcs = ["testFeatures.cpp"],
    data = [
        "lactotroph.cl",
        "//clode/cpp:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//clode/cpp",
    ] + select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["@com_github_khronosgroup_opencl//:opencl_headers"],
    }),
)

cc_test(
    name = "testFeatLarge",
    srcs = ["testFeatures_morevars.cpp"],
    data = [
        ":lactotroph",
        "//clode/cpp:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//clode/cpp",
    ] + select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["@com_github_khronosgroup_opencl//:opencl_headers"],
    }),
)

cc_test(
    name = "listOpenCL",
    srcs = ["listOpenCL.cpp"],
    data = [
        "//clode/cpp:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//clode/cpp",
    ] + select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": [],
        "@platforms//os:windows": ["@com_github_khronosgroup_opencl//:opencl_headers"],
    }),
)