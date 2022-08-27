cc_library(
    name = "opencl_resources",
    srcs = ["OpenCLResource.cpp"],
    hdrs = ["OpenCLResource.hpp"],
    deps = [
        "//src/OpenCL:cl2",
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
)

cc_library(
    name = "clode_api",
    hdrs = [
        "CLODE.hpp",
        "clODE_random.cl",
        "clODE_struct_defs.cl",
        "clODE_utilities.cl",
    ],
    deps = [
        ":cl_definitions",
        "//src/OpenCL:cl2",
    ],
)

cc_library(
    name = "clode",
    srcs = ["CLODE.cpp"],
    deps = [
        ":clode_api",
        ":opencl_resources",
        "//src/OpenCL:cl2",
    ],
)
