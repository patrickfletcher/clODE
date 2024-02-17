import clode


def test_print_open_cl_features():
    input_file: str = "test/van_der_pol_oscillator.cl"

    t_span = (0.0, 1000.0)

    features = clode.FeatureSimulator(
        src_file=input_file,
        variables={"x": 1.0, "y": 1.0},
        parameters={"mu": 1.0},
        num_noise=0,
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.dormand_prince,
        t_span=t_span,
    )

    features.print_devices()


def test_print_open_cl_trajectory():
    input_file: str = "test/van_der_pol_oscillator.cl"

    t_span = (0.0, 1000.0)

    trajectory = clode.TrajectorySimulator(
        src_file=input_file,
        variables={"x": 1.0, "y": 1.0},
        parameters={"mu": 1.0},
        num_noise=0,
        stepper=clode.Stepper.dormand_prince,
        t_span=t_span,
    )

    trajectory.print_devices()
