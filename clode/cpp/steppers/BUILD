package(default_visibility = ["//visibility:public"])

filegroup(
    name = "cl",
    srcs = [
        "adaptive_bs23.clh",
        "adaptive_dp45.clh",
        "adaptive_explicit_step.clh",
        "adaptive_he12.clh",
        "adaptive_rkck.clh",
        "fixed_explicit_Euler.clh",
        "fixed_explicit_RK4.clh",
        "fixed_explicit_Trapezoidal.clh",
        "fixed_explicit_step.clh",
        "fixed_explicit_stochastic_Euler.clh",
    ],
)

cc_library(
    name = "steppers",
    hdrs = [
        "adaptive_bs23.clh",
        "adaptive_dp45.clh",
        "adaptive_explicit_step.clh",
        "adaptive_he12.clh",
        "adaptive_rkck.clh",
        "fixed_explicit_Euler.clh",
        "fixed_explicit_RK4.clh",
        "fixed_explicit_Trapezoidal.clh",
        "fixed_explicit_step.clh",
        "fixed_explicit_stochastic_Euler.clh",
    ],
    deps = [
        "//clode/cpp:cl_definitions",
    ],
)
