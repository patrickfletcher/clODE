load("//bazel:repositories.bzl", "clode_dependencies")

clode_dependencies()

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "rules_python",
    sha256 = "b593d13bb43c94ce94b483c2858e53a9b811f6f10e1e0eedc61073bd90e58d9c",
    strip_prefix = "rules_python-0.12.0",
    url = "https://github.com/bazelbuild/rules_python/archive/refs/tags/0.12.0.tar.gz",
)

http_archive(
    name = "build_bazel_rules_apple",
    sha256 = "90e3b5e8ff942be134e64a83499974203ea64797fd620eddeb71b3a8e1bff681",
    url = "https://github.com/bazelbuild/rules_apple/releases/download/1.1.2/rules_apple.1.1.2.tar.gz",
)

load("@rules_python//python:repositories.bzl", "python_register_toolchains")

python_register_toolchains(
    name = "python3_9",
    # Available versions are listed in @rules_python//python:versions.bzl.
    # We recommend using the same version your team is already standardized on.
    python_version = "3.9",
)
#
#load("@python3_9//:defs.bzl", "interpreter")
#load("@rules_python//python:pip.bzl", "pip_parse")

#pip_parse(
#    ...
#    python_interpreter_target = interpreter,
#    ...
#)

#load(
#    "@build_bazel_rules_apple//apple:repositories.bzl",
#    "apple_rules_dependencies",
#)
#
#apple_rules_dependencies()
#
#load(
#    "@build_bazel_rules_swift//swift:repositories.bzl",
#    "swift_rules_dependencies",
#)
#
#swift_rules_dependencies()
#
#load(
#    "@build_bazel_rules_swift//swift:extras.bzl",
#    "swift_rules_extra_dependencies",
#)
#
#swift_rules_extra_dependencies()
#
#load(
#    "@build_bazel_apple_support//lib:repositories.bzl",
#    "apple_support_dependencies",
#)
#
#apple_support_dependencies()

load("//bazel/python:python_native.bzl", "clode_python_deps")

clode_python_deps()
