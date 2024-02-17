from typing import List

from clode import OpenCLConverter
import math

def get_rhs(t: float,
        x: List[float],
        p: List[float],
        dx: List[float],
        aux: List[float],
        w: List[float]) -> None:
    v: float = x[0]
    n: float = x[1]
    c: float = x[2]

    gca: float = p[0]
    gkca: float = p[1]
    kpmca: float = p[2]
    gk: float = 3500.0

    vca: float = 25.0
    vk: float = -75.0
    cm: float = 5300.0
    alpha: float = 4.5e-6
    fcyt: float = 0.01
    kd: float = 0.4
    vm: float = -20.0
    sm: float = 12.0
    vn: float = -16.0
    sn: float = 5.0
    taun: float = 20.0

    minf: float = 1.0 / (1.0 + math.exp((vm - v) / sm))
    ninf: float = 1.0 / (1.0 + math.exp((vn - v) / sn))
    omega: float = c ** 2 / (c ** 2 + kd ** 2)

    ica: float = gca * minf * (v - vca)
    ik: float = gk * n * (v - vk)
    ikca: float = gkca * omega * (v - vk)

    dx[0] = -(ica + ik + ikca) / cm
    dx[1] = (ninf - n) / taun
    dx[2] = fcyt * (-alpha * ica - kpmca * c)


converter = OpenCLConverter()
print(converter.convert_to_opencl(get_rhs))
print("Foo")
# import clode
# from clode import OpenCLConverter, OpenCLExp
#
#
# def helper_function(x: float, y: float) -> float:
#     return (x - y) / OpenCLExp(y)
#
#
# def get_rhs(
#     t: float,
#     variables: list[float],
#     parameters: list[float],
#     derivatives: list[float],
#     aux: list[float],
#     wiener: list[float],
# ) -> None:
#     x: float = variables[0]
#     y: float = variables[1]
#
#     derivatives[0] = helper_function(x, y)
#     derivatives[1] = -x
#
#
# converter = OpenCLConverter()
# _ = converter.convert_to_opencl(helper_function)
# opencl = converter.convert_to_opencl(get_rhs)
#
# simulator = clode.TrajectorySimulator(
#     rhs_equation=get_rhs,
#     variables={"x": 0.0, "y": 2.0},
#     parameters={"a": 1.0},
#     supplementary_equations=[helper_function],
#     t_span=(0, 10),
#     dt=0.01,
# )
#
# trajs = simulator.trajectory()
#
#
# print(opencl)
