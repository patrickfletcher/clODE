from pathlib import Path

import pytest

from clode._backends.rhs import create_rhs_source, load_rhs_source
from clode._pyopencl import (
    KernelKind,
    KernelRegistry,
    Precision,
    ProblemShape,
    RhsValidationError,
    SourceBuilder,
    UnsupportedObserverError,
    UnsupportedStepperError,
)
from test.core_numerics.helpers import model_path


KERNEL_ROOT = Path(__file__).resolve().parents[1] / "clode" / "cpp"


def test_kernel_registry_exposes_current_cpp_defines_and_entrypoints() -> None:
    registry = KernelRegistry(KERNEL_ROOT)

    assert registry.get_stepper_define("rk4") == "EXPLICIT_RK4"
    assert registry.get_observer_define("basicall") == "USE_OBSERVER_BASIC_ALLVAR"
    assert tuple(path.name for path in registry.get_entrypoint_paths(KernelKind.FEATURES)) == (
        "transient.cl",
        "initializeObserver.cl",
        "features.cl",
    )
    assert registry.get_kernel_names(KernelKind.TRAJECTORY) == (
        "transient",
        "trajectory",
    )


def test_kernel_registry_rejects_unknown_stepper_and_observer() -> None:
    registry = KernelRegistry(KERNEL_ROOT)

    with pytest.raises(UnsupportedStepperError):
        registry.get_stepper_define("bogus")

    with pytest.raises(UnsupportedObserverError):
        registry.get_observer_define("bogus")


def test_source_builder_matches_transient_source_assembly_and_build_options() -> None:
    builder = SourceBuilder(KERNEL_ROOT)
    rhs = load_rhs_source(model_path("stable_linear.cl"))

    bundle = builder.build(
        kernel_kind=KernelKind.TRANSIENT,
        precision=Precision.SINGLE,
        stepper_name="rk4",
        problem_shape=ProblemShape(n_var=1, n_par=1, n_aux=0, n_wiener=0),
        rhs=rhs,
    )

    assert bundle.build_key.kernel_kind is KernelKind.TRANSIENT
    assert bundle.build_key.rhs_digest == rhs.digest
    assert bundle.build_options == (
        "-DCLODE_SINGLE_PRECISION",
        "-DEXPLICIT_RK4",
        "-DN_PAR=1",
        "-DN_VAR=1",
        "-DN_AUX=0",
        "-DN_WIENER=0",
        f"-I{KERNEL_ROOT}",
    )
    assert bundle.kernel_names == ("transient",)
    assert "__kernel void transient" in bundle.source_text
    assert bundle.source_text.endswith(rhs.text)


def test_source_builder_features_include_observer_options_and_build_key_changes() -> None:
    builder = SourceBuilder(KERNEL_ROOT)
    rhs_a = create_rhs_source("rhs_a.cl", "void getRHS() {}\n")
    rhs_b = create_rhs_source("rhs_b.cl", "void getRHS() {\n    return;\n}\n")

    first = builder.build(
        kernel_kind=KernelKind.FEATURES,
        precision=Precision.DOUBLE,
        stepper_name="dopri5",
        problem_shape=ProblemShape(n_var=2, n_par=1, n_aux=1, n_wiener=0),
        rhs=rhs_a,
        observer_name="basicall",
        n_store_events=7,
        debug_build=True,
    )
    second = builder.build(
        kernel_kind=KernelKind.FEATURES,
        precision=Precision.DOUBLE,
        stepper_name="dopri5",
        problem_shape=ProblemShape(n_var=2, n_par=1, n_aux=1, n_wiener=0),
        rhs=rhs_b,
        observer_name="basicall",
        n_store_events=7,
        debug_build=True,
    )

    assert first.kernel_names == ("transient", "initializeObserver", "features")
    assert "__kernel void initializeObserver" in first.source_text
    assert "__kernel void features" in first.source_text
    assert "-DUSE_OBSERVER_BASIC_ALLVAR" in first.build_options
    assert "-DN_STORE_EVENTS=7" in first.build_options
    assert "-g" in first.build_options
    assert "-cl-opt-disable" in first.build_options
    assert first.build_key.rhs_digest != second.build_key.rhs_digest
    assert first.build_key.kernel_tree_digest == second.build_key.kernel_tree_digest


def test_source_builder_validates_rhs_and_feature_configuration() -> None:
    builder = SourceBuilder(KERNEL_ROOT)

    with pytest.raises(RhsValidationError):
        builder.build(
            kernel_kind=KernelKind.TRANSIENT,
            precision=Precision.SINGLE,
            stepper_name="rk4",
            problem_shape=ProblemShape(n_var=1, n_par=0, n_aux=0, n_wiener=0),
            rhs=create_rhs_source("invalid_rhs.cl", "void notRHS() {}\n"),
        )

    with pytest.raises(UnsupportedObserverError):
        builder.build(
            kernel_kind=KernelKind.FEATURES,
            precision=Precision.SINGLE,
            stepper_name="rk4",
            problem_shape=ProblemShape(n_var=1, n_par=0, n_aux=0, n_wiener=0),
            rhs=create_rhs_source("valid_rhs.cl", "void getRHS() {}\n"),
        )

    with pytest.raises(ValueError, match="observer_name"):
        builder.build(
            kernel_kind=KernelKind.TRANSIENT,
            precision=Precision.SINGLE,
            stepper_name="rk4",
            problem_shape=ProblemShape(n_var=1, n_par=0, n_aux=0, n_wiener=0),
            rhs=create_rhs_source("valid_rhs.cl", "void getRHS() {}\n"),
            observer_name="basic",
        )