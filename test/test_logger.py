import clode
from test.core_numerics.helpers import device_kwargs_for_tests


def test_print_open_cl(capfd):
    input_file: str = "test/van_der_pol_oscillator.cl"

    trajectory = clode.TrajectorySimulator(
        src_file=input_file,
        variables={"x": 0.0, "y": 1.0},
        parameters={"mu": 1.0},
        num_noise=0,
        stepper=clode.Stepper.dormand_prince,
        **device_kwargs_for_tests(platform_id=0, device_id=0),
    )

    clode.set_log_level(clode.LogLevel.trace)
    assert clode.get_log_level() == clode.LogLevel.trace
    trajectory.print_devices()
    captured = capfd.readouterr()
    assert "OpenCL" in captured.out
    assert captured.err == ""
    clode.set_log_level(clode.LogLevel.off)
    assert clode.get_log_level() == clode.LogLevel.off
    trajectory.print_devices()
    captured = capfd.readouterr()
    assert captured.out == ""
    assert captured.err == ""
