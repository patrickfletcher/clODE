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
    name = "test_pituitary",
    size = "large",
    timeout = "short",
    srcs = ["test_pituitary.py"],
    data = [":test.cl"],
    deps = [
        "//pyclode:clode_init",
        "//src:py_clode",
    ],
)

py_test(
    name = "test_vdp",
    srcs = ["test_vdp.py"],
    data = [":van_der_pol_oscillator.cl"],
    deps = [
        "//pyclode:clode_init",
        "//src:py_clode",
    ],
)

py_test(
    name = "test_vdp_long",
    size = "large",
    timeout = "short",
    srcs = [
        "test_vdp.py",
        "test_vdp_long.py",
    ],
    data = [":van_der_pol_oscillator.cl"],
    deps = [
        "//pyclode:clode_init",
        "//src:py_clode",
    ],
)
