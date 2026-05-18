from __future__ import annotations

import numpy as np
import pytest

import clode
from test.core_numerics.helpers import STABLE_LINEAR_PARAMETERS, STABLE_LINEAR_VARIABLES, model_path


def affine_rhs(
    t: float,
    variables: list[float],
    parameters: list[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
) -> None:
    rate: float = parameters[0]
    offset: float = parameters[1]
    derivatives[0] = rate * variables[0] + offset


def coupled_rhs(
    t: float,
    variables: list[float],
    parameters: list[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
) -> None:
    a: float = parameters[0]
    b: float = parameters[1]
    x: float = variables[0]
    y: float = variables[1]

    derivatives[0] = a * x + b
    derivatives[1] = x - y
    aux[0] = x + y


def test_initial_value_problem_construction_preserves_defaults_and_metadata() -> None:
    ivp = clode.InitialValueProblem(
        variables={"x": 1.5, "y": -2.0},
        parameters={"a": 0.5, "b": 1.25},
        aux=["sum"],
        rhs_equation=coupled_rhs,
    )

    assert ivp.variable_names == ["x", "y"]
    assert ivp.parameter_names == ["a", "b"]
    assert ivp.aux_names == ["sum"]
    assert ivp.num_variables == 2
    assert ivp.num_parameters == 2
    assert ivp.num_aux == 1
    assert ivp.num_noise == 0
    assert ivp.ensemble_size == 1
    assert ivp.ensemble_shape == (1,)
    assert ivp.python_rhs_equation is coupled_rhs
    assert "void getRHS" in ivp.rhs_source.text
    np.testing.assert_allclose(ivp.default_initial_state, [1.5, -2.0])
    np.testing.assert_allclose(ivp.default_parameters, [0.5, 1.25])


def test_initial_value_problem_batch_shaping_keeps_shape_and_default_fallback() -> None:
    ivp = clode.InitialValueProblem(
        variables={"x": 1.0, "y": 2.0},
        parameters={"a": 3.0, "b": 4.0},
        aux=["sum"],
        rhs_equation=coupled_rhs,
    )

    ivp.set_ensemble(
        variables={"x": np.array([[1.0, 2.0], [3.0, 4.0]])},
        parameters={"b": np.array([[10.0, 20.0], [30.0, 40.0]])},
    )

    assert ivp.ensemble_size == 4
    assert ivp.ensemble_shape == (2, 2)
    np.testing.assert_allclose(
        ivp.get_initial_state(),
        [[1.0, 2.0], [2.0, 2.0], [3.0, 2.0], [4.0, 2.0]],
    )
    np.testing.assert_allclose(
        ivp.get_parameter_values(),
        [[3.0, 10.0], [3.0, 20.0], [3.0, 30.0], [3.0, 40.0]],
    )

    ivp.set_repeat_ensemble(3)

    assert ivp.ensemble_size == 3
    assert ivp.ensemble_shape == (3, 1)
    np.testing.assert_allclose(ivp.get_initial_state(), np.tile([[1.0, 2.0]], (3, 1)))
    np.testing.assert_allclose(ivp.get_parameter_values(), np.tile([[3.0, 4.0]], (3, 1)))


def test_python_backed_ivp_call_supports_default_and_scipy_args_order() -> None:
    ivp = clode.InitialValueProblem(
        variables={"x": 0.0},
        parameters={"rate": -0.5, "offset": 1.0},
        rhs_equation=affine_rhs,
    )

    np.testing.assert_allclose(ivp(0.0, np.array([0.0])), [1.0])
    np.testing.assert_allclose(ivp(0.0, np.array([2.0]), -0.25, 0.5), [0.0])
    np.testing.assert_allclose(
        ivp.evaluate_rhs(0.0, np.array([2.0]), parameters={"offset": 0.0}),
        [-1.0],
    )

    with pytest.raises(ValueError, match="Expected 2 positional parameter overrides"):
        ivp(0.0, np.array([2.0]), -0.25)


def test_non_python_backed_ivp_is_not_callable() -> None:
    ivp = clode.InitialValueProblem(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
    )

    with pytest.raises(ValueError, match="Python RHS"):
        ivp(0.0, np.array([1.0, 2.0]))


def test_python_backed_ivp_works_with_solve_ivp_and_args() -> None:
    scipy_integrate = pytest.importorskip("scipy.integrate")

    ivp = clode.InitialValueProblem(
        variables={"x": 0.0},
        parameters={"rate": -0.5, "offset": 1.0},
        rhs_equation=affine_rhs,
    )

    solution = scipy_integrate.solve_ivp(
        ivp,
        (0.0, 1.0),
        [0.0],
        args=(-0.5, 1.0),
        rtol=1e-10,
        atol=1e-12,
    )

    expected = 2.0 * (1.0 - np.exp(-0.5))
    assert solution.success
    np.testing.assert_allclose(solution.y[:, -1], [expected], rtol=1e-9, atol=1e-11)