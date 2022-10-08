#load("@com_github_grpc_grpc//third_party/py:python_configure.bzl", "python_configure")
load(":python.bzl", "python_configure")

def clode_python_deps():
    python_configure(name = "local_config_python")

    native.bind(
        name = "python_headers",
        actual = "@local_config_python//:python_headers",
    )
