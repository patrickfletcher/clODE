import clode


def test_print_open_cl(capfd):
    input_file: str = "test/van_der_pol_oscillator.cl"

    tspan = (0.0, 1000.0)

    trajectory = clode.CLODETrajectory(
        src_file=input_file,
        variable_names=["x", "y"],
        parameter_names=["mu"],
        num_noise=0,
        stepper=clode.Stepper.dormand_prince,
        tspan=tspan,
        device_id=0,
        platform_id=0,
    )

    clode.set_log_level(clode.log_level.trace)
    assert clode.get_log_level() == clode.log_level.trace
    trajectory.print_devices()
    captured = capfd.readouterr()
    assert "OpenCL" in captured.out
    assert captured.err == ""
    clode.set_log_level(clode.log_level.off)
    assert clode.get_log_level() == clode.log_level.off
    trajectory.print_devices()
    captured = capfd.readouterr()
    assert captured.out == ""
    assert captured.err == ""
