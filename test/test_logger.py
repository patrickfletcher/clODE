from dataclasses import replace
import io
import logging
import os
from pathlib import Path
import warnings

import clode
import pytest

from clode._opencl import BuildError, KernelKind, OpenCLRuntime, Precision, ProblemShape, SourceBuilder
from clode.problem._core import load_rhs_source
from clode.runtime import _clode_root_dir
from test.core_numerics.helpers import device_kwargs_for_tests, model_path


pyopencl = pytest.importorskip("pyopencl")
KERNEL_ROOT = Path(_clode_root_dir)


def _reset_logging_state() -> None:
    logging.captureWarnings(False)

    clode_logger = logging.getLogger("clode")
    for handler in list(clode_logger.handlers):
        clode_logger.removeHandler(handler)
        handler.close()
    clode_logger.addHandler(logging.NullHandler())
    clode_logger.setLevel(logging.WARNING)
    clode_logger.propagate = True

    warnings_logger = logging.getLogger("py.warnings")
    for handler in list(warnings_logger.handlers):
        warnings_logger.removeHandler(handler)
        handler.close()
    warnings_logger.setLevel(logging.NOTSET)
    warnings_logger.propagate = True


def _transient_source_bundle():
    builder = SourceBuilder(KERNEL_ROOT)
    return builder.build(
        kernel_kind=KernelKind.TRANSIENT,
        precision=Precision.SINGLE,
        stepper_name="rk4",
        problem_shape=ProblemShape(n_var=2, n_par=2, n_aux=0, n_wiener=0),
        rhs=load_rhs_source(model_path("stable_linear.cl")),
    )


def test_print_devices_reports_visible_opencl_devices(capfd):
    trajectory = clode.TrajectorySimulator(
        src_file="test/van_der_pol_oscillator.cl",
        variables={"x": 0.0, "y": 1.0},
        parameters={"mu": 1.0},
        num_noise=0,
        stepper=clode.Stepper.dormand_prince,
        **device_kwargs_for_tests(platform_id=0, device_id=0),
    )

    trajectory.print_devices()
    captured = capfd.readouterr()
    assert "OpenCL" in captured.out
    assert captured.err == ""


def test_configure_logging_routes_clode_records_and_pyopencl_warnings(monkeypatch):
    stream = io.StringIO()

    try:
        monkeypatch.delenv("PYOPENCL_COMPILER_OUTPUT", raising=False)

        clode.configure_logging(
            level="INFO",
            stream=stream,
            capture_warnings=True,
            compiler_output=True,
        )

        logger = clode.get_logger("runtime.audit")
        logger.debug("hidden debug")
        logger.info("visible info")
        warnings.warn("compiler output available", pyopencl.CompilerWarning)

        output = stream.getvalue()
        assert "visible info" in output
        assert "hidden debug" not in output
        assert "compiler output available" in output
        assert os.environ["PYOPENCL_COMPILER_OUTPUT"] == "1"
    finally:
        _reset_logging_state()


def test_build_error_mentions_pyopencl_compiler_output_hint(monkeypatch):
    monkeypatch.delenv("PYOPENCL_COMPILER_OUTPUT", raising=False)

    runtime = OpenCLRuntime.create(**device_kwargs_for_tests(platform_id=0, device_id=0))
    source_bundle = _transient_source_bundle()
    invalid_bundle = replace(
        source_bundle,
        source_text=source_bundle.source_text + "\nthis is not valid opencl source\n",
    )

    with pytest.raises(BuildError) as exc_info:
        runtime.program_cache.get_or_build(runtime, invalid_bundle)

    assert "PYOPENCL_COMPILER_OUTPUT=1" in exc_info.value.build_log
