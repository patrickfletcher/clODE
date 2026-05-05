from pathlib import Path

import clode

from clode._backends.rhs import compute_rhs_digest
from test.core_numerics.helpers import (
    STABLE_LINEAR_PARAMETERS,
    STABLE_LINEAR_VARIABLES,
    device_kwargs_for_tests,
    model_path,
)


FIXED_DT = 0.05
FIXED_MAX_STEPS = 64


def test_file_rhs_source_is_prepared_with_text_and_digest() -> None:
    source_path = model_path("stable_linear.cl")
    expected_text = Path(source_path).read_text(encoding="utf-8")

    simulator = clode.Simulator(
        src_file=source_path,
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator._rhs_source.origin_label == source_path
    assert simulator._rhs_source.text == expected_text
    assert simulator._rhs_source.digest == compute_rhs_digest(expected_text)
    assert expected_text in simulator.get_program_string()


def test_python_rhs_source_is_prepared_with_text_and_digest() -> None:
    def get_rhs(
        t: float,
        state: list[float],
        parameters: list[float],
        derivatives: list[float],
        aux: list[float],
        noise: list[float],
    ) -> None:
        derivatives[0] = -parameters[0] * state[0]

    simulator = clode.Simulator(
        rhs_equation=get_rhs,
        variables={"x": 1.0},
        parameters={"k": 1.0},
        aux=[],
        num_noise=0,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    assert simulator._rhs_source.origin_label == "clode_rhs.cl"
    assert "void getRHS" in simulator._rhs_source.text
    assert simulator._rhs_source.digest == compute_rhs_digest(simulator._rhs_source.text)
    assert simulator._rhs_source.text in simulator.get_program_string()


def test_rhs_digest_changes_when_source_text_changes() -> None:
    first = compute_rhs_digest("void getRHS() { return; }\n")
    second = compute_rhs_digest("void getRHS() {\n    return;\n}\n")

    assert first != second
