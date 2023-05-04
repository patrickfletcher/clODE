import clode


def test_print_open_cl_features():
    input_file: str = "test/van_der_pol_oscillator.cl"

    tspan = (0.0, 1000.0)

    features = clode.CLODEFeatures(
        src_file=input_file,
        variable_names=["x", "y"],
        parameter_names=["mu"],
        num_noise=0,
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.dormand_prince,
        tspan=tspan,
    )

    features.print_devices()


def test_print_open_cl_trajectory():
    input_file: str = "test/van_der_pol_oscillator.cl"

    tspan = (0.0, 1000.0)

    trajectory = clode.CLODETrajectory(
        src_file=input_file,
        variable_names=["x", "y"],
        parameter_names=["mu"],
        num_noise=0,
        stepper=clode.Stepper.dormand_prince,
        tspan=tspan,
    )

    trajectory.print_devices()
