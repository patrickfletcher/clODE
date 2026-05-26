from __future__ import annotations

from collections.abc import Iterable

import pytest


MARKERS_BY_PREFIX: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("test/core_numerics/", ("core_numerics", "numerics", "release_gate", "requires_opencl")),
    ("test/kernel_components/", ("kernel_components", "requires_opencl")),
)

MARKERS_BY_FILE: dict[str, tuple[str, ...]] = {
    "test/test_function_converter.py": ("frontend", "smoke"),
    "test/test_initial_value_problem.py": ("frontend", "smoke"),
    "test/test_initial_value_problem_runtime.py": ("runtime_api", "release_gate", "requires_opencl"),
    "test/test_opencl_builtins.py": ("frontend", "release_gate", "requires_opencl"),
    "test/test_xpp_parser.py": ("frontend", "release_gate", "requires_opencl"),
    "test/test_simulation_contracts.py": ("runtime_api", "release_gate", "requires_opencl"),
    "test/test_problem_rhs_source.py": ("runtime_api", "release_gate", "requires_opencl"),
    "test/test_opencl_runtime.py": ("runtime_api", "release_gate", "requires_opencl"),
    "test/test_runtime.py": ("runtime_api", "release_gate", "requires_opencl"),
    "test/test_logger.py": ("runtime_api", "release_gate", "requires_opencl"),
    "test/test_ornl_thompson_a1.py": ("numerics", "requires_opencl"),
    "test/test_vdp.py": ("numerics", "requires_opencl"),
    "test/test_aux_values.py": ("numerics", "requires_opencl"),
    "test/test_features.py": ("numerics", "requires_opencl"),
    "test/test_opencl_models.py": ("opencl_internal", "smoke"),
    "test/test_opencl_source_builder.py": ("opencl_internal", "smoke"),
    "test/test_opencl_buffers.py": ("opencl_internal", "requires_opencl"),
    "test/test_opencl_structs.py": ("opencl_internal", "requires_opencl"),
}


def _markers_for_path(path: str) -> Iterable[str]:
    for prefix, markers in MARKERS_BY_PREFIX:
        if path.startswith(prefix):
            yield from markers
    yield from MARKERS_BY_FILE.get(path, ())


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    for item in items:
        path = item.path.as_posix()
        for marker_name in _markers_for_path(path):
            item.add_marker(getattr(pytest.mark, marker_name))
