from __future__ import annotations

import subprocess
import sys


PRIMARY_BUNDLES: dict[str, tuple[str, ...]] = {
    "smoke": (
        "test/test_function_converter.py",
        "test/test_pyopencl_models.py",
        "test/test_pyopencl_source_builder.py",
    ),
    "frontend": (
        "test/test_function_converter.py",
        "test/test_opencl_builtins.py",
        "test/test_xpp_parser.py",
    ),
    "runtime_api": (
        "test/test_backend_contracts.py",
        "test/test_backend_rhs_source.py",
        "test/test_pyopencl_runtime.py",
        "test/test_runtime.py",
        "test/test_logger.py",
    ),
    "numerics": (
        "test/core_numerics",
        "test/test_ornl_thompson_a1.py",
        "test/test_vdp.py",
        "test/test_aux_values.py",
        "test/test_features.py",
    ),
    "pyopencl_internal": (
        "test/test_pyopencl_models.py",
        "test/test_pyopencl_source_builder.py",
        "test/test_pyopencl_buffers.py",
        "test/test_pyopencl_structs.py",
    ),
    "legacy_cpp_comparison": (
        "test/test_pyopencl_transient_backend.py",
        "test/test_pyopencl_trajectory_backend.py",
        "test/test_pyopencl_feature_backend.py",
    ),
}

ALIASES: dict[str, tuple[str, ...]] = {
    "release": ("frontend", "runtime_api", "numerics"),
    "opencl": ("release",),
    "extended": ("pyopencl_internal", "legacy_cpp_comparison"),
}

DESCRIPTIONS = {
    "smoke": "Cross-platform packaging and driver-independent smoke checks.",
    "frontend": "Frontend-facing conversion and API-surface tests.",
    "runtime_api": "Public runtime, backend selection, and simulator contract tests.",
    "numerics": "Numerical regression tests, including the core OpenCL reference suite.",
    "pyopencl_internal": "Focused tests for internal PyOpenCL support layers.",
    "legacy_cpp_comparison": "Comparison tests that still rely on the legacy C++ wrapper backend.",
    "release": "The OpenCL-backed release gate: frontend, runtime/API, and numerical regressions.",
    "opencl": "Compatibility alias for the release gate bundle.",
    "extended": "Additional internal and legacy-comparison checks for manual verification.",
}


def _resolve_bundle(bundle_name: str) -> tuple[str, ...]:
    if bundle_name in PRIMARY_BUNDLES:
        return PRIMARY_BUNDLES[bundle_name]
    if bundle_name in ALIASES:
        ordered_paths: list[str] = []
        for alias_target in ALIASES[bundle_name]:
            for candidate in _resolve_bundle(alias_target):
                if candidate not in ordered_paths:
                    ordered_paths.append(candidate)
        return tuple(ordered_paths)
    raise KeyError(bundle_name)


def _print_usage() -> None:
    print("Usage: python tools/run_test_bundle.py <bundle> [-- <pytest args>]")
    print("")
    print("Available bundles:")
    for name in sorted((*PRIMARY_BUNDLES.keys(), *ALIASES.keys())):
        print(f"  {name}: {DESCRIPTIONS[name]}")
        for path in _resolve_bundle(name):
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
        test_paths = _resolve_bundle(bundle_name)
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
