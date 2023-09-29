from math import exp

import numpy as np

import clode


def ornl_thompson_a1_exact(t: float):
    y1 = 4 * (t + 1 / 8 * exp(-8 * t) - 1 / 8)
    y2 = 4 * (1 - exp(-8 * t))
    return y1, y2


def test_ornl_thompson_a1():
    num_simulations = 1

    m = 1 / 4
    w = 8
    k = 2
    H = 10

    tspan = (0.0, H / 2.0)

    integrator = clode.CLODETrajectory(
        src_file="test/ornl_thompson_a1.cl",
        variable_names=["y1", "y2"],
        parameter_names=["m", "w", "k", "H"],
        aux=["g1"],
        num_noise=0,
        dt=0.001,
        dtmax=0.001,
        stepper=clode.Stepper.rk4,
        tspan=tspan,
    )

    parameters = [m, w, k, H]

    x0 = np.tile([0, 0], (num_simulations, 1))
    pars_v = np.tile(parameters, (num_simulations, 1))

    integrator.initialize(x0, pars_v)

    integrator.trajectory()

    time_steps = integrator.get_time_steps()
    output_trajectory = integrator.get_trajectory()[0]

    for tt, (y1, y2) in zip(time_steps, output_trajectory):
        expected_y1, expected_y2 = ornl_thompson_a1_exact(tt)
        np.testing.assert_approx_equal(y1, expected_y1, significant=5)
        np.testing.assert_approx_equal(y2, expected_y2, significant=5)


# if using 'bazel test ...'
if __name__ == "__main__":
    test_ornl_thompson_a1()
