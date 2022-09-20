cc_library(
    name = "matlab_sink",
    hdrs = ["MatlabSink.hpp"],
    deps = [
        "@com_github_gabime_spdlog//:spdlog",
    ],
)

cc_library(
    name = "python_sink",
    hdrs = ["PythonSink.hpp"],
    visibility = ["//visibility:public"],
    deps = [
        "@com_github_gabime_spdlog//:spdlog",
        "@com_github_pybind_pybind11//:pybind11",
    ],
)
