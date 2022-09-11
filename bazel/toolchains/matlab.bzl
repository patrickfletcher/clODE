load("/bazel/toolchains/mac/matlab_utils.bzl", "configure_osx_mex_toolchain")



mex_library = rule(
    implementation = _mex_library_impl,
    attrs = {
        "deps": attr.label_list(),
        ...
    },
)


def _mex_library_impl(ctx):
