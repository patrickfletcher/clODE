filegroup(
    name = "lactotroph",
    srcs = ["test.cl"],
    visibility = ["//visibility:public"],
)

py_test(
    name = "test_lactotroph",
    srcs = ["test_lactotroph.py"],
    deps = [
        "//clode/python:clode_all_sources",
    ],
)

py_test(
    name = "test_lactotroph_more_vars",
    size = "large",
    timeout = "short",
    srcs = ["test_lactotroph_more_vars.py"],
    data = [":test.cl"],
    deps = [
        "//clode/python:clode_all_sources",
    ],
)

py_test(
    name = "test_lactotroph_trajectory",
    srcs = ["test_lactotroph_trajectory.py"],
    data = [
        ":test.cl",
        "//samples:lactotroph",
    ],
    deps = [
        "//clode/python:clode_all_sources",
    ],
)

py_test(
    name = "test_pituitary",
    size = "large",
    timeout = "short",
    srcs = ["test_pituitary.py"],
    data = [":test.cl"],
    deps = [
        "//clode/python:clode_all_sources",
    ],
)

py_test(
    name = "test_vdp",
    srcs = ["test_vdp.py"],
    data = [":van_der_pol_oscillator.cl"],
    deps = [
        "//clode/python:clode_all_sources",
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
        "//clode/python:clode_all_sources",
    ],
)

py_test(
    name = "test_ornl_thompson_a1",
    #size = "large",
    #timeout = "short",
    srcs = [
        "test_ornl_thompson_a1.py",
    ],
    data = [":ornl_thompson_a1.cl"],
    deps = [
        "//clode/python:clode_all_sources",
    ],
)
