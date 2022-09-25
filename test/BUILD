#py_library(
#    name = "clode_pylib",
#    data = ["//src:clode.so"],
#    imports = ["."],
#    visibility = ["//visibility:private"],
#)

py_test(
    name = "test_clode",
    srcs = ["test_clode.py"],
    deps = ["//src:py_clode"],
)

py_test(
    name = "pituitary",
    srcs = ["pituitary.py"],
    data = [":test.cl"],
    deps = ["//src:py_clode"],
)
