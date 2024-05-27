from __future__ import annotations

import sys
from math import log, pi

import numpy as np
import pytest

import clode


def cuberoot(x):
    return x ** (1 / 3.0)


def approximate_vdp_period(mu):
    # https://www.johndcook.com/blog/2019/12/26/van-der-pol-period/
    # https://math.stackexchange.com/questions/1564464/how-to-find-the-period-of-periodic-solutions-of-the-van-der-pol-equation
    if mu < 0:
        period = 0
    elif 0 <= mu < 2:
        period = 2 * pi * (1 + mu**2 / 16)
    else:
        period = min(
            2 * pi * (1 + mu**2 / 16),
            (3 - 2 * log(2)) * mu + 3 * 2.2338 / cuberoot(mu),
        )
    return period


def getRHS(
    t: float,
    var: list[float],
    par: list[float],
    derivatives: list[float],
    aux: list[float],
    wiener: list[float],
) -> None:
    mu: float = par[0]
    x: float = var[0]
    y: float = var[1]

    dx: float = y
    dy: float = mu * (1 - x * x) * y - x

    derivatives[0] = dx
    derivatives[1] = dy


def vdp_dormand_prince(
    end: int,
    input_file: str | None = None,
    input_eq: clode.OpenCLRhsEquation | None = None,
):
    if input_file is None and input_eq is None:
        input_file = "test/van_der_pol_oscillator.cl"
    t_span = (0.0, 1000.0)

    integrator = clode.FeatureSimulator(
        src_file=input_file,
        rhs_equation=input_eq,
        variables={"x": 1.0, "y": 1.0},
        parameters={"mu": 1.0},
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.dormand_prince,
        t_span=t_span,
        max_store=20000,
        max_steps=20000,
    )

    parameters = {
        "mu": [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
        + list(range(5, end))
    }

    integrator.set_ensemble(parameters=parameters)

    # integrator.transient()
    integrator.features()
    observer_output = integrator.get_observer_results()

    periods = observer_output.get_var_max("period")

    for index, mu in enumerate(parameters["mu"]):
        period = periods[index]
        expected_period = approximate_vdp_period(mu)
        rtol = 0.01
        atol = 1
        assert np.isclose(period, expected_period, rtol=rtol, atol=atol), (
            f"Period {period} not close to expected {expected_period}" + f"for mu {mu}"
        )


def test_vdp_dormand_prince():
    vdp_dormand_prince(end=7)


def test_vdp_dormand_prince_python_rhs():
    vdp_dormand_prince(end=7, input_eq=getRHS)


# if using 'bazel test ...'
if __name__ == "__main__":
    # print(clode)
    # sys.exit(pytest.main(sys.argv[1:]))
    # clode.set_log_level(clode.LogLevel.debug)
    vdp_dormand_prince(end=100)
