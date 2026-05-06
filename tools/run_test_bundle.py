from __future__ import annotations

import subprocess
import sys


BUNDLES: dict[str, tuple[str, ...]] = {
    "smoke": (
        "test/test_function_converter.py",
        "test/test_pyopencl_models.py",
        "test/test_pyopencl_source_builder.py",
    ),
    "opencl": (
        "test/core_numerics",
        "test/test_backend_contracts.py",
        "test/test_backend_rhs_source.py",
        "test/test_pyopencl_runtime.py",
        "test/test_runtime.py",
        "test/test_logger.py",
    ),
    "extended": (
        "test/test_pyopencl_buffers.py",
        "test/test_pyopencl_structs.py",
        "test/test_pyopencl_transient_backend.py",
        "test/test_pyopencl_trajectory_backend.py",
        "test/test_pyopencl_feature_backend.py",
        "test/test_aux_values.py",
        "test/test_features.py",
        "test/test_trajectory.py",
        "test/test_opencl_builtins.py",
        "test/test_vdp.py",
        "test/test_xpp_parser.py",
    ),
    "long": (
        "test/test_vdp_long.py",
    ),
}

DESCRIPTIONS = {
    "smoke": "Cross-platform package and pure-Python regression checks.",
    "opencl": "Linux OpenCL runtime gate used by default CI.",
    "extended": "Additional runtime and integration checks for manual verification.",
    "long": "Opt-in long-running tests.",
}


def _print_usage() -> None:
    print("Usage: python tools/run_test_bundle.py <bundle> [-- <pytest args>]")
    print("")
    print("Available bundles:")
    for name in sorted(BUNDLES):
        print(f"  {name}: {DESCRIPTIONS[name]}")
        for path in BUNDLES[name]:
            print(f"    - {path}")


def main(argv: list[str]) -> int:
    if not argv or argv[0] in {"-h", "--help"}:
        _print_usage()
        return 0

    if argv[0] == "--list":
        _print_usage()
        return 0

    bundle_name = argv[0]
    try:
        test_paths = BUNDLES[bundle_name]
    except KeyError:
        print(f"Unknown bundle: {bundle_name}", file=sys.stderr)
        print("", file=sys.stderr)
        _print_usage()
        return 2

    pytest_args = argv[1:]
    if pytest_args and pytest_args[0] == "--":
        pytest_args = pytest_args[1:]

    command = [sys.executable, "-m", "pytest", "-q", *test_paths, *pytest_args]
    print("Running test bundle:", bundle_name)
    return subprocess.call(command)


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
