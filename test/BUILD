#py_library(
#    name = "clode_pylib",
#    data = ["//src:clode.so"],
#    imports = ["."],
#    visibility = ["//visibility:private"],
#)

py_test(
    name = "test_lactotroph",
    srcs = ["test_lactotroph.py"],
    deps = ["//src:py_clode"],
)

py_test(
    name = "test_lactotroph_more_vars",
    srcs = ["test_lactotroph_more_vars.py"],
    data = [":test.cl"],
    deps = ["//src:py_clode"],
)

py_test(
    name = "pituitary",
    srcs = ["pituitary.py"],
    data = [":test.cl"],
    deps = ["//src:py_clode"],
)
