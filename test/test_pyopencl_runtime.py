from dataclasses import replace
from pathlib import Path

import pytest

pyopencl = pytest.importorskip("pyopencl")

from clode._pyopencl import BuildError, KernelKind, OpenCLRuntime, Precision, ProblemShape, SourceBuilder
from test.core_numerics.helpers import TEST_DEVICE_ID, TEST_PLATFORM_ID, model_path


KERNEL_ROOT = Path(__file__).resolve().parents[1] / "clode" / "cpp"


def _explicit_runtime_kwargs() -> dict[str, int]:
    return {
        "platform_id": 0 if TEST_PLATFORM_ID is None else TEST_PLATFORM_ID,
        "device_id": 0 if TEST_DEVICE_ID is None else TEST_DEVICE_ID,
    }


def _transient_source_bundle():
    from clode._backends.rhs import load_rhs_source

    builder = SourceBuilder(KERNEL_ROOT)
    return builder.build(
        kernel_kind=KernelKind.TRANSIENT,
        precision=Precision.SINGLE,
        stepper_name="rk4",
        problem_shape=ProblemShape(n_var=2, n_par=2, n_aux=0, n_wiener=0),
        rhs=load_rhs_source(model_path("stable_linear.cl")),
    )


def test_opencl_runtime_selects_explicit_device_and_reports_capabilities() -> None:
    runtime = OpenCLRuntime.create(**_explicit_runtime_kwargs())

    assert runtime.platform_id == _explicit_runtime_kwargs()["platform_id"]
    assert runtime.device_id == _explicit_runtime_kwargs()["device_id"]
    assert runtime.get_max_memory_alloc_size() > 0
    assert "OpenCL" in runtime.get_device_cl_version()


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