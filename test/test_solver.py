from typing import List

import numpy as np
import pytest

import clode
from clode import ObserverParams, SolverParams

"""
Test the Simulator class API
- set/get that move data between host/device, not invalidating CL program
- appropriate manipulation of initial condition and parameter arrays to/from expected device layout
"""


def get_rhs(
    t: float,
    vars: List[float],
    p: List[float],
    dy: List[float],
    aux: List[float],
    w: List[float],
) -> None:
    a: float = p[0]
    b: float = p[1]
    c: float = p[2]
    x: float = vars[0]
    y: float = vars[1]
    dy[0] = y
    dy[1] = -a * x
    aux[0] = b * c * x


variables = {"x": -1.0, "y": 0.0}
parameters = {"a": 1.0, "b": 1.0, "c": -0.5}
aux = ["aux_x"]
num_noise = 0

op = clode.ObserverParams()

simulator = clode.Simulator(
    rhs_equation=get_rhs,
    parameters=parameters,
    variables=variables,
    aux=aux,
    num_noise=num_noise,
)


@pytest.mark.skip(reason="not ready yet")
def test_set_get_nobuild():

    # print(simulator.get_tspan())
    # simulator.set_tspan((0, 15))
    # print(simulator.get_tspan())
    # simulator.shift_tspan()
    # print(simulator.get_tspan())

    # print(simulator.get_solver_parameters())
    # simulator.set_solver_parameters(dt=1e-3, dtmax=10, abstol=1e-5, reltol=1e-4, max_steps=10000, max_store=1000, nout=10)
    # print(simulator.get_solver_parameters())
    # sp = clode.SolverParams(dt=1e-2, abstol=2e-5,)
    # simulator.set_solver_parameters(sp)
    # print(simulator.get_solver_parameters())

    # simulator.print_devices()
    # simulator.print_status()
    pass


@pytest.mark.skip(reason="not ready yet")
def test_set_ensemble():

    # TODO: fixtures
    init_vars = variables
    # init_vars = variables
    # init_vars = np.zeros((3,2))
    # pars = parameters
    # pars={"a": [[0.5, 1.0],[2.0, 3.0]], "b": 1.1, "c": np.linspace(0, 5.,4)}
    pars = {
        "a": [[0.5, 1.0], [1.5, 2], [2.5, 3.0]],
        "b": 1.1,
        "c": [[0, 1], [2, 3], [4, 5]],
    }
    simulator.set_ensemble(variables=init_vars, parameters=pars)
    print(simulator._device_initial_state)
    print(simulator._device_parameters)
    print(simulator._ensemble_size, simulator._ensemble_shape)

    pass


@pytest.mark.skip(reason="not ready yet")
def test_set_repeat_ensemble():
    # simulator.set_repeat_ensemble(3)
    # print(simulator.get_initial_state())
    # print(simulator._device_parameters)
    pass


def run_tests():
    test_set_get_nobuild()
    test_set_ensemble()
    test_set_repeat_ensemble()


if __name__ == "__main__":
    run_tests()
