load(":clode_http_archive.bzl", "clode_http_archive")
load(":external_deps.bzl", "load_repository_locations")
load(":repository_locations.bzl", "REPOSITORY_LOCATIONS_SPEC")

REPOSITORY_LOCATIONS = load_repository_locations(REPOSITORY_LOCATIONS_SPEC)

def external_http_archive(name, **kwargs):
    clode_http_archive(
        name,
        locations = REPOSITORY_LOCATIONS,
        **kwargs
    )

def clode_dependencies(skip_targets = []):
    _com_github_fmtlib_fmt()
    _com_github_gabime_spdlog()
    _org_python_ftp()
    _com_github_pybind11_pybind()
    _com_github_khronosgroup_opencl_windows()

def _com_github_fmtlib_fmt():
    external_http_archive(
        name = "com_github_fmtlib_fmt",
        build_file = "//bazel/external:fmtlib.BUILD",
    )
    native.bind(
        name = "fmtlib",
        actual = "@com_github_fmtlib_fmt//:fmtlib",
    )

def _com_github_gabime_spdlog():
    external_http_archive(
        name = "com_github_gabime_spdlog",
        build_file = "//bazel/external:spdlog.BUILD",
    )
    native.bind(
        name = "spdlog",
        actual = "@com_github_gabime_spdlog//:spdlog",
    )

def _org_python_ftp():
    external_http_archive(
        name = "org_python_ftp",
        build_file = "//bazel/external:python.BUILD",
    )
    native.bind(
        name = "spdlog",
        actual = "@org_python_ftp//:python",
    )

def _com_github_pybind11_pybind():
    external_http_archive(
        name = "com_github_pybind_pybind11",
        build_file = "//bazel/external:pybind11.BUILD",
    )
    native.bind(
        name = "pybind11",
        actual = "@com_github_pybind_pybind11//:pybind11",
    )

def _com_github_khronosgroup_opencl_windows():
    external_http_archive(
        name = "com_github_khronosgroup_opencl_windows",
        build_file = "//bazel/external:opencl_windows.BUILD",
    )
    native.bind(
        name = "opencl_windows",
        actual = "@com_github_khronosgroup_opencl_windows//:opencl_windows",
    )

#def _com_github_grpc_grpc():
#    external_http_archive(
#        name = "com_github_grpc_grpc",
#        build_file = "//bazel/external:grpc.BUILD",
#    )
#    native.bind(
#        name = "grpc",
#        actual = "@com_github_grpc_grpc//:grpc",
#    )
