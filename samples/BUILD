cc_binary(
    name = "testTrans",
    srcs = ["testTransient.cpp"],
    data = [
        "lactotroph.cl",
        "//src:cl_definitions",
        "//src:clode_cl_api",
        "//src:opencl",
    ],
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_binary(
    name = "testTraj",
    srcs = ["testTrajectory.cpp"],
    data = [
        "lactotroph.cl",
        "//src:cl_definitions",
        "//src:clode_cl_api",
        "//src:opencl",
    ],
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_binary(
    name = "testFeat",
    srcs = ["testFeatures.cpp"],
    data = [
        "lactotroph.cl",
        "//src:opencl",
    ],
    linkstatic = True,
    deps = [
        "//src",
    ],
)

cc_binary(
    name = "listOpenCL",
    srcs = ["listOpenCL.cpp"],
    data = [
        "//src:opencl",
    ],
    linkstatic = True,
    deps = [
        "//src",
    ],
)
