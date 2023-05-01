import clode


def test_print_open_cl(capfd):
    clode.set_log_level(clode.log_level.trace)
    assert clode.get_log_level() == clode.log_level.trace
    clode.print_devices()
    captured = capfd.readouterr()
    assert "OpenCL" in captured.out
    assert captured.err == ""
    clode.set_log_level(clode.log_level.off)
    assert clode.get_log_level() == clode.log_level.off
    clode.print_devices()
    captured = capfd.readouterr()
    assert captured.out == ""
    assert captured.err == ""
