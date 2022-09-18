package(default_visibility = ["//visibility:public"])

load(":matlab_toolchain_config.bzl", "cc_toolchain_config")

cc_toolchain_config(name = "mex_toolchain_config")

filegroup(name = "empty")

filegroup(
    name = "all",
    srcs = [
        "mex_osx.sh",
    ],
)

cc_toolchain(
    name = "mex_toolchain",
    all_files = ":all",
    compiler_files = ":all",
    dwp_files = ":all",
    linker_files = ":all",
    objcopy_files = ":all",
    strip_files = ":all",
    supports_param_files = 0,
    toolchain_config = ":mex_toolchain_config",
    toolchain_identifier = "mex-toolchain",
)

cc_toolchain_suite(
    name = "matlab_mex",
    toolchains = {
        "mex": ":mex_toolchain",
    },
)
