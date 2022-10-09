#py_library(
#    name = "clode_pylib",
#    data = ["//src:clode.so"],
#    imports = ["."],
#    visibility = ["//visibility:private"],
#)

filegroup(
    name = "lactotroph",
    srcs = ["test.cl"],
    visibility = ["//visibility:public"],
)

py_test(
    name = "test_lactotroph",
    srcs = ["test_lactotroph.py"],
    deps = ["//src:py_clode"],
)

py_test(
    name = "test_lactotroph_more_vars",
    size = "large",
    timeout = "short",
    srcs = ["test_lactotroph_more_vars.py"],
    data = [":test.cl"],
    deps = ["//src:py_clode"],
)

py_test(
    name = "test_lactotroph_trajectory",
    srcs = ["test_lactotroph_trajectory.py"],
    data = [
        ":test.cl",
        "//samples:lactotroph",
    ],
    deps = ["//src:py_clode"],
)

py_test(
    name = "pituitary",
    size = "large",
    timeout = "short",
    srcs = ["pituitary.py"],
    data = [":test.cl"],
    deps = ["//src:py_clode"],
)
