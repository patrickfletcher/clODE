import pytest
from typing import List
import clode
from clode import SolverParams, ObserverParams

'''
Test the Simulator class API
- setters that cause cl_initialized=False (requiring CL program rebuild)
- set/get that move data between host/device, not invalidating CL program
- appropriate manipulation of initial condition and parameter arrays to/from expected device layout
'''

def get_rhs(
    t: float,
    vars: List[float],
    p: List[float],
    dy: List[float],
    aux: List[float],
    w: List[float],
) -> None:
    k: float = p[0]
    x: float = vars[0]
    y: float = vars[1]
    dy[0] = y
    dy[1] = -k * x
    aux[0] = x

variables = {"x": -1.0, "y": 0.0}
parameters = {"k": 1.0}
aux = ["aux_x"]
num_noise = 0

op = clode.ObserverParams()

simulator = clode.Simulator(rhs_equation=get_rhs, parameters=parameters, variables=variables, aux=aux, num_noise=num_noise)

def test_set_get_nobuild():

    print(simulator.get_tspan())
    simulator.set_tspan((0, 15))
    print(simulator.get_tspan())
    simulator.shift_tspan()
    print(simulator.get_tspan())
    print(simulator.is_initialized)

    print(simulator.get_solver_parameters())
    simulator.set_solver_parameters(dt=1e-3, dtmax=10, abstol=1e-5, reltol=1e-4, max_steps=10000, max_store=1000, nout=10)
    print(simulator.get_solver_parameters())
    sp = clode.SolverParams(dt=1e-2, abstol=2e-5,)
    simulator.set_solver_parameters(sp)
    print(simulator.get_solver_parameters())
    print(simulator.is_initialized)
    
    print(simulator.get_initial_state())
    print(simulator.get_final_state())

    simulator.print_devices()
    simulator.print_status()

def run_tests():
    test_set_get_nobuild()

if __name__=="__main__":
    run_tests()