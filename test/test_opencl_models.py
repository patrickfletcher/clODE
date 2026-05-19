from dataclasses import replace

import numpy as np
import pytest

from clode._opencl import (
    BuildError,
    BuildKey,
    DoublePrecisionNotSupportedError,
    KernelKind,
    OPENCL_BACKEND_VERSION,
    Precision,
    ProblemShape,
    ProgramBundle,
    RhsValidationError,
    SourceBundle,
    UnsupportedStepperError,
)
from clode.observers._definitions import resolve_observer_spec
from clode.observers import ObserverParams, get_observer_feature_names, is_two_pass_observer
from clode.problem._core import ProblemInfo, create_rhs_source


def _make_build_key(**overrides: object) -> BuildKey:
    rhs = create_rhs_source("model.cl", "void getRHS() {}\n")
    build_key_kwargs = {
        "backend_version": OPENCL_BACKEND_VERSION,
        "kernel_kind": KernelKind.TRANSIENT,
        "precision": Precision.SINGLE,
        "stepper_name": "rk4",
        "observer_define": None,
        "problem_shape": ProblemShape(n_var=2, n_par=1, n_aux=1, n_wiener=0),
        "n_store_events": 0,
        "rhs_digest": rhs.digest,
        "kernel_tree_digest": "kernel-tree-digest",
        "debug_build": False,
    }
    build_key_kwargs.update(overrides)
    return BuildKey(**build_key_kwargs)


def test_problem_shape_from_problem_info_matches_solver_dimensions() -> None:
    problem_info = ProblemInfo(
        "stable_linear.cl",
        ["x", "y"],
        ["k"],
        ["aux0"],
        2,
    )

    assert ProblemShape.from_problem_info(problem_info) == ProblemShape(
        n_var=2,
        n_par=1,
        n_aux=1,
        n_wiener=2,
    )


def test_build_key_is_hashable_and_validates_observer_event_contract() -> None:
    build_key = _make_build_key(precision=Precision.from_single_precision(True))

    assert {build_key: "cached"}[build_key] == "cached"

    with pytest.raises(ValueError, match="observer_define"):
        _make_build_key(n_store_events=1)


def test_program_bundle_validates_build_key_and_kernel_handles() -> None:
    build_key = _make_build_key()
    source_bundle = SourceBundle(
        build_key=build_key,
        source_text="__kernel void transient(void) {}\n",
        build_options=("-DCLODE_SINGLE_PRECISION",),
        kernel_names=("transient",),
    )
    kernel_handle = object()

    bundle = ProgramBundle(
        build_key=build_key,
        source_bundle=source_bundle,
        program=object(),
        kernels={"transient": kernel_handle},
    )

    assert bundle.kernels["transient"] is kernel_handle

    with pytest.raises(ValueError, match="source bundle"):
        ProgramBundle(
            build_key=replace(build_key, debug_build=True),
            source_bundle=source_bundle,
            program=object(),
            kernels={"transient": object()},
        )

    with pytest.raises(ValueError, match="missing kernel handles"):
        ProgramBundle(
            build_key=build_key,
            source_bundle=source_bundle,
            program=object(),
            kernels={},
        )


def test_opencl_errors_preserve_diagnostic_context() -> None:
    build_error = BuildError(
        "build failed",
        source_text="kernel source",
        build_options=("-DTEST=1", "-I/tmp/kernels"),
        build_log="compiler output",
    )

    assert build_error.source_text == "kernel source"
    assert build_error.formatted_build_options == "-DTEST=1 -I/tmp/kernels"
    assert build_error.build_log == "compiler output"

    rhs_error = RhsValidationError("missing getRHS", origin_label="generated_rhs.cl")
    stepper_error = UnsupportedStepperError("bogus-stepper")
    precision_error = DoublePrecisionNotSupportedError("Mock GPU")

    assert rhs_error.origin_label == "generated_rhs.cl"
    assert "generated_rhs.cl" in str(rhs_error)
    assert stepper_error.stepper_name == "bogus-stepper"
    assert "bogus-stepper" in str(stepper_error)
    assert precision_error.device_name == "Mock GPU"
    assert "Mock GPU" in str(precision_error)


def test_observer_catalog_helpers_expose_feature_names_and_two_pass_flags() -> None:
    problem_info = ProblemInfo(
        "stable_linear_aux.cl",
        ["x", "y"],
        ["k"],
        ["aux0"],
        0,
    )

    feature_names = get_observer_feature_names(
        problem_info,
        "basic",
        ObserverParams(f_var_ix=1),
    )

    assert feature_names == (
        "max y",
        "min y",
        "mean y",
        "max dy/dt",
        "min dy/dt",
        "step count",
    )
    assert is_two_pass_observer("nhood2") is True
    assert is_two_pass_observer("basic") is False


def test_resolved_observer_spec_splits_persistent_and_event_layout() -> None:
    problem_info = ProblemInfo(
        "stable_linear_aux.cl",
        ["x", "y"],
        ["k"],
        ["aux0"],
        0,
    )

    resolved_spec = resolve_observer_spec(
        problem_info,
        "localmax",
        ObserverParams(max_event_timestamps=2),
        real_dtype=np.dtype(np.float32),
    )

    assert resolved_spec.observer_name == "localmax"
    assert resolved_spec.build_define == "USE_OBSERVER_LOCAL_MAX"
    assert resolved_spec.uses_two_pass is False
    assert tuple(field[0] for field in resolved_spec.layout.event_output_fields) == (
        "tMaxList",
        "xMaxList",
        "tMinList",
        "xMinList",
    )
    assert "tMaxList" not in tuple(
        field[0] for field in resolved_spec.layout.persistent_fields
    )
    assert resolved_spec.observer_data_struct_name == (
        "clode_observer_data_localmax_float_v2_a1_e2"
    )