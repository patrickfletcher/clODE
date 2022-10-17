filegroup(
    name = "lactotroph",
    srcs = ["lactotroph.cl"],
    visibility = ["//visibility:public"],
)

cc_test(
    name = "testTrans",
    srcs = ["testTransient.cpp"],
    data = [
        "lactotroph.cl",
        "//src:cl_definitions",
        "//src:clode_cl_api",
        "//src:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_test(
    name = "testTraj",
    srcs = ["testTrajectory.cpp"],
    data = [
        "lactotroph",
        "//src:cl_definitions",
        "//src:clode_cl_api",
        "//src:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_test(
    name = "testFeat",
    srcs = ["testFeatures.cpp"],
    data = [
        "lactotroph.cl",
        "//src:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_test(
    name = "testFeatLarge",
    srcs = ["testFeatures_morevars.cpp"],
    data = [
        "//src:opencl",
        "//test:lactotroph",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_test(
    name = "listOpenCL",
    srcs = ["listOpenCL.cpp"],
    data = [
        "//src:opencl",
    ],
    linkopts = select({
        "@platforms//os:osx": [],
        "@platforms//os:linux": ["-lOpenCL"],
        "@platforms//os:windows": [],
    }),
    linkstatic = True,
    deps = [
        "//src",
    ],
)
