import pytest
from math import cos, exp, pi, sqrt
from typing import List

import numpy as np

import clode


def sine_curve(
    t: float,
    x_: List[float],
    p_: List[float],
    dx_: List[float],
    aux_: List[float],
    w_: List[float],
) -> None:
    x: float = x_[0]
    dilation: float = p_[0]
    dx: float = cos(t * dilation)
    dx_[0] = dx
    aux_[0] = x + 1.0
    aux_[1] = 1.0
    aux_[2] = -2.0


@pytest.mark.parametrize("observer",
                         [clode.Observer.basic_all_variables,
                          clode.Observer.local_max,
                          clode.Observer.neighbourhood_1,
                          clode.Observer.neighbourhood_2,
                          clode.Observer.threshold_2]
)
def test_sine_curve_timestamps(observer: clode.Observer):
    "Test that the active timestamps of a sine curve are correct"

    # Define the parameters
    parameters = {
        "dilation": 1.0,
    }

    # Define the initial state
    variables = {
        "x": 0.0,
    }

    # Activate at t=pi/4, deactivate at t=3pi/2
    feature_simulator = clode.FeatureSimulator(
        rhs_equation=sine_curve,
        variables=variables,
        parameters=parameters,
        aux=["xp1", "pos", "neg"],
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.rk4,
        dtmax=0.1,
        dt=0.1,
        t_span=(0.0, 400 * pi),
        event_var="x",
        feature_var="x",
    )

    # Run the simulation
    output = feature_simulator.features()

    assert output is not None

    # Check that the min, max and mean aux values are correct

    xp1_mean = output.get_var_mean("xp1")
    xp1_min = output.get_var_min("xp1")
    xp1_max = output.get_var_max("xp1")

    assert xp1_mean == pytest.approx(1.0, abs=1e-2)
    assert xp1_min == pytest.approx(0.0, abs=1e-2)
    assert xp1_max == pytest.approx(2.0, abs=1e-2)

    pos_mean = output.get_var_mean("pos")
    pos_min = output.get_var_min("pos")
    pos_max = output.get_var_max("pos")

    assert pos_mean == 1.0
    assert pos_min == 1.0
    assert pos_max == 1.0

    neg_mean = output.get_var_mean("neg")
    neg_min = output.get_var_min("neg")
    neg_max = output.get_var_max("neg")

    # assert neg_mean == -1.2
    assert neg_min == -2.0
    assert neg_max == -2.0

    pass