from math import exp

import pyclode
import pyclode.stepper


def ornl_thompson_a1_exact(t: float):
    y1 = 4 * (t * 1 / 8 * (exp(-8 * t) - 1))
    y2 = 4 * (1 - exp(-8 * t))
    return y1, y2


def test_ornl_thompson_a1():

    m = 1 / 4
    w = 8
    k = 2
    H = 10

    tspan = (0.0, H / 2)

    integrator = pyclode.CLODETrajectory(
        src_file="test/van_der_pol_oscillator.cl",
        variable_names=["y1", "y2"],
        parameter_names=["m", "w", "k", "H"],
        aux_names=["g1"],
        num_noise=0,
        observer=pyclode.Observer.threshold_2,
        stepper=pyclode.stepper.Stepper.dormand_prince,
        tspan=tspan,
    )

    parameters = [m, w, k, H]

    x0 = np.tile([1, 1], (len(parameters), 1))
    pars_v = np.tile(parameters, (len(parameters), 1))

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
        assert np.isclose(period, expected_period, rtol=rtol, atol=1), \
            f"Period {period} not close to expected {expected_period}" + \
            f"for mu {mu}"
