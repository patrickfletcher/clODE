#py_library(
#    name = "clode_pylib",
#    data = ["//src:clode.so"],
#    imports = ["."],
#    visibility = ["//visibility:private"],
#)

py_test(
    name = "test_clode",
    srcs = ["test_clode.py"],
    data = [
        "//samples:lactotroph",
        "//src:opencl",
    ],
    deps = ["//src:py_clode"],
)
