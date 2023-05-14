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


def vdp_dormand_prince(end: int, input_file: str = "test/van_der_pol_oscillator.cl"):
    tspan = (0.0, 1000.0)

    integrator = clode.CLODEFeatures(
        src_file=input_file,
        variable_names=["x", "y"],
        parameter_names=["mu"],
        num_noise=0,
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.dormand_prince,
        tspan=tspan,
        single_precision=False,
    )

    parameters = [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] + list(
        range(5, end)
    )

    x0 = np.tile([1, 1], (len(parameters), 1))

    pars_v = np.array([[par] for par in parameters])
    integrator.initialize(x0, pars_v)

    integrator.transient()
    integrator.features()
    observer_output = integrator.get_observer_results()

    periods = observer_output.get_var_max("period")
    for index, mu in enumerate(parameters):
        period = periods[index, 0]
        expected_period = approximate_vdp_period(mu)
        rtol = 0.01
        atol = 0.3
        assert np.isclose(period, expected_period, rtol=rtol, atol=1), (
            f"Period {period} not close to expected {expected_period}" + f"for mu {mu}"
        )


def test_vdp_dormand_prince():
    vdp_dormand_prince(end=7)


# if using 'bazel test ...'
if __name__ == "__main__":
    print(clode)
    sys.exit(pytest.main(sys.argv[1:]))
