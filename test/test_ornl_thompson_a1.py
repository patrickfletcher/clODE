from math import exp

import numpy as np

import clode

# test type: numerical accuracy (compare simulation to exact solution)


def ornl_thompson_a1_exact(t: float):
    y1 = 4.0 * (t + 1.0 / 8.0 * exp(-8.0 * t) - 1.0 / 8.0)
    y2 = 4.0 * (1.0 - exp(-8.0 * t))
    return y1, y2


def test_ornl_thompson_a1():

    m = 1 / 4
    w = 8.0
    k = 2.0
    H = 10.0

    t_span = (0.0, H / 2.0)

    parameters = {"m": m, "w": w, "k": k, "H": H}

    integrator = clode.TrajectorySimulator(
        src_file="test/ornl_thompson_a1.cl",
        variables={"y1": 0.0, "y2": 0.0},
        parameters=parameters,
        aux=["g1"],
        num_noise=0,
        dt=0.001,
        dtmax=0.001,
        stepper=clode.Stepper.rk4,
        t_span=t_span,
        max_store=20000,
        max_steps=20000,
    )

    trajectory = integrator.trajectory()

    for tt, y1, y2 in zip(trajectory.t, trajectory.x["y1"], trajectory.x["y2"]):
        expected_y1, expected_y2 = ornl_thompson_a1_exact(tt)
        np.testing.assert_approx_equal(y1, expected_y1, significant=5)
        np.testing.assert_approx_equal(y2, expected_y2, significant=5)


# if using 'bazel test ...'
if __name__ == "__main__":
    test_ornl_thompson_a1()
