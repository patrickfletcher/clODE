from dataclasses import replace
import json
import os
from pathlib import Path
import subprocess
import sys
import textwrap
from types import SimpleNamespace

import pytest

import clode

pyopencl = pytest.importorskip("pyopencl")

from clode._opencl import BuildError, KernelKind, OpenCLRuntime, Precision, ProblemShape, SourceBuilder
from clode._opencl.runtime import _selection_sort_key
from clode.runtime.selection import RuntimeSelection
from clode.problem._core import load_rhs_source
from clode.runtime import _clode_root_dir
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID, model_path


KERNEL_ROOT = Path(_clode_root_dir)


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def _transient_source_bundle():
    builder = SourceBuilder(KERNEL_ROOT)
    return builder.build(
        kernel_kind=KernelKind.TRANSIENT,
        precision=Precision.SINGLE,
        stepper_name="rk4",
        problem_shape=ProblemShape(n_var=2, n_par=2, n_aux=0, n_wiener=0),
        rhs=load_rhs_source(model_path("stable_linear.cl")),
    )


def _fake_platform(name: str, vendor: str) -> SimpleNamespace:
    return SimpleNamespace(name=name, vendor=vendor)


def _fake_device(device_type: int, *, vendor: str, name: str) -> SimpleNamespace:
    return SimpleNamespace(type=device_type, vendor=vendor, name=name)


def test_opencl_runtime_selects_explicit_device_and_reports_capabilities() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())

    assert runtime.platform_id == _explicit_runtime_kwargs()["platform_id"]
    assert runtime.device_id == _explicit_runtime_kwargs()["device_id"]
    assert runtime.get_max_memory_alloc_size() > 0
    assert "OpenCL" in runtime.get_device_cl_version()


def test_opencl_runtime_default_selection_prefers_gpu_when_available() -> None:
    platforms = pyopencl.get_platforms()
    has_gpu = any(
        int(device.type) & int(clode.CLDeviceType.DEVICE_TYPE_GPU)
        for platform in platforms
        for device in platform.get_devices()
    )
    if not has_gpu:
        pytest.skip("No GPU runtime is visible in this environment")

    runtime = OpenCLRuntime.create()

    assert int(runtime.device.type) & int(clode.CLDeviceType.DEVICE_TYPE_GPU)


def test_selection_sort_key_prefers_gpu_over_cpu() -> None:
    gpu_key = _selection_sort_key(
        1,
        0,
        _fake_platform("NVIDIA CUDA", "NVIDIA Corporation"),
        _fake_device(
            int(clode.CLDeviceType.DEVICE_TYPE_GPU),
            vendor="NVIDIA Corporation",
            name="NVIDIA GeForce RTX 5080",
        ),
    )
    cpu_key = _selection_sort_key(
        0,
        0,
        _fake_platform("Portable Computing Language", "The pocl project"),
        _fake_device(
            int(clode.CLDeviceType.DEVICE_TYPE_CPU),
            vendor="The pocl project",
            name="pthread-cpu",
        ),
    )

    assert gpu_key < cpu_key


def test_selection_sort_key_prefers_non_pocl_cpu_over_pocl_cpu() -> None:
    regular_cpu_key = _selection_sort_key(
        0,
        0,
        _fake_platform("Intel(R) OpenCL", "Intel(R) Corporation"),
        _fake_device(
            int(clode.CLDeviceType.DEVICE_TYPE_CPU),
            vendor="Intel(R) Corporation",
            name="Intel(R) Xeon",
        ),
    )
    pocl_cpu_key = _selection_sort_key(
        1,
        0,
        _fake_platform("Portable Computing Language", "The pocl project"),
        _fake_device(
            int(clode.CLDeviceType.DEVICE_TYPE_CPU),
            vendor="The pocl project",
            name="pthread-cpu",
        ),
    )

    assert regular_cpu_key < pocl_cpu_key


def test_opencl_runtime_from_selection_accepts_explicit_device_id() -> None:
    runtime = OpenCLRuntime.from_selection(
        RuntimeSelection(
            device_type=None,
            vendor=None,
            platform_id=_explicit_runtime_kwargs()["platform_id"],
            device_id=_explicit_runtime_kwargs()["device_id"],
        )
    )

    assert runtime.platform_id == _explicit_runtime_kwargs()["platform_id"]
    assert runtime.device_id == _explicit_runtime_kwargs()["device_id"]


def test_simulator_exposes_concrete_selected_runtime() -> None:
    simulator = clode.TrajectorySimulator(
        src_file="test/van_der_pol_oscillator.cl",
        variables={"x": 0.0, "y": 1.0},
        parameters={"mu": 1.0},
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 10.0),
        **_explicit_runtime_kwargs(),
    )

    assert simulator.platform_id == _explicit_runtime_kwargs()["platform_id"]
    assert simulator.device_id == _explicit_runtime_kwargs()["device_id"]
    assert simulator.runtime_selection == RuntimeSelection(
        device_type=None,
        vendor=None,
        platform_id=_explicit_runtime_kwargs()["platform_id"],
        device_id=_explicit_runtime_kwargs()["device_id"],
    )
    assert "platform_id=" in simulator.runtime_description


def test_opencl_runtime_from_selection_requires_explicit_device_id() -> None:
    with pytest.raises(ValueError, match="Must specify device_id"):
        OpenCLRuntime.from_selection(
            RuntimeSelection(
                device_type=None,
                vendor=None,
                platform_id=_explicit_runtime_kwargs()["platform_id"],
                device_id=None,
            )
        )


def test_program_cache_is_runtime_scoped_and_reuses_build_key_hits() -> None:
    source_bundle = _transient_source_bundle()
    first_runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    second_runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())

    first_bundle = first_runtime.program_cache.get_or_build(first_runtime, source_bundle)
    cached_bundle = first_runtime.program_cache.get_or_build(first_runtime, source_bundle)
    second_bundle = second_runtime.program_cache.get_or_build(second_runtime, source_bundle)

    assert first_bundle is cached_bundle
    assert len(first_runtime.program_cache) == 1
    assert len(second_runtime.program_cache) == 1
    assert first_bundle is not second_bundle
    assert set(first_bundle.kernels) == {"transient"}


def test_program_cache_surfaces_build_failures_with_source_and_options() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())
    source_bundle = _transient_source_bundle()
    invalid_bundle = replace(
        source_bundle,
        source_text=source_bundle.source_text + "\nthis is not valid opencl source\n",
    )

    with pytest.raises(BuildError) as exc_info:
        runtime.program_cache.get_or_build(runtime, invalid_bundle)

    error = exc_info.value
    assert error.source_text == invalid_bundle.source_text
    assert error.build_options == invalid_bundle.build_options
    assert error.build_log.strip() != ""


def test_default_public_path_does_not_import_cpp_wrapper() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = textwrap.dedent(
        f"""
        import builtins
        import json
        import os
        import sys

        attempts = []
        real_import = builtins.__import__

        def guarded_import(name, globals=None, locals=None, fromlist=(), level=0):
            if name == \"clode.cpp\" or name.startswith(\"clode.cpp.\"):
                attempts.append({{\"name\": name, \"fromlist\": list(fromlist or ())}})
                raise ModuleNotFoundError(name)
            return real_import(name, globals, locals, fromlist, level)

        builtins.__import__ = guarded_import

        import clode

        platforms = clode.query_opencl()
        simulator = clode.Simulator(
            src_file=\"test/van_der_pol_oscillator.cl\",
            variables={{\"x\": 0.0, \"y\": 1.0}},
            parameters={{\"mu\": 1.0}},
            num_noise=0,
            stepper=clode.Stepper.rk4,
            platform_id={_explicit_runtime_kwargs()["platform_id"]},
            device_id={_explicit_runtime_kwargs()["device_id"]},
        )

        print(
            json.dumps(
                {{
                    \"attempts\": attempts,
                    \"platform_count\": len(platforms),
                    \"program_string_len\": len(simulator.get_program_string()),
                    \"wrapper_loaded\": \"clode.cpp.clode_cpp_wrapper\" in sys.modules,
                }}
            )
        )
        """
    )
    completed = subprocess.run(
        [sys.executable, "-c", script],
        check=True,
        capture_output=True,
        cwd=repo_root,
        env=os.environ.copy(),
        text=True,
    )
    payload = json.loads(completed.stdout.strip().splitlines()[-1])

    assert payload["platform_count"] > 0
    assert payload["program_string_len"] > 0
    assert payload["wrapper_loaded"] is False
    assert payload["attempts"] == []