from typing import List

import numpy as np

import clode
from clode import (
    TrajectorySimulator,
    acos,
    acosh,
    acospi,
    asin,
    asinh,
    asinpi,
    atan,
    atan2,
    atan2pi,
    atanh,
    atanpi,
    cbrt,
    ceil,
    copysign,
    cos,
    cosh,
    cospi,
    erf,
    erfc,
    exp,
    exp2,
    exp10,
    expm1,
    fabs,
    fdim,
    floor,
    fmod,
    gamma,
    heaviside,
    hypot,
    ilogb,
    ldexp,
    lgamma,
    log,
    log1p,
    log2,
    log10,
    nextafter,
    pow,
    pown,
    powr,
    remainder,
    rint,
    rootn,
    rsqrt,
    sin,
    sinh,
    sinpi,
    sqrt,
    tan,
    tanh,
    tanpi,
    trunc,
)


def test_opencl_builtins() -> None:
    def rhs(
        t: float,
        x: List[float],
        p: List[float],
        dx: List[float],
        aux: List[float],
        w: List[float],
    ) -> None:
        #  Test every OpenCL builtin function and place the result in aux

        p0: int = int(p[0])
        p1: int = int(p[1])
        x0: float = x[0] * t
        x1: float = x[1] * t
        x2: float = x[2] * t

        aux[0] = acos(x0)
        aux[1] = acosh(x0 + 1)
        aux[2] = acospi(x0)
        aux[3] = asin(x0)
        aux[4] = asinh(x0 + 1)
        aux[5] = asinpi(x0)
        aux[6] = atan(x0)
        aux[7] = atan2(x0, x1)
        aux[8] = atan2pi(x0, x1)
        aux[9] = atanh(x0 / 100.0)
        aux[10] = atanpi(x0)
        aux[11] = cbrt(x0)
        aux[12] = ceil(x0)
        aux[13] = copysign(x0, x1)
        aux[14] = cos(x0)
        aux[15] = cosh(x0)
        aux[16] = cospi(x0)
        aux[17] = erf(x0)
        aux[18] = erfc(x0)
        aux[19] = exp(x0)
        aux[20] = exp2(x0)
        aux[21] = exp10(x0)
        aux[22] = expm1(x0)
        aux[23] = fabs(x0)
        aux[24] = fdim(x0, x1)
        aux[25] = floor(x1)
        aux[26] = 0.0  # fma(x0, x1, x2)
        aux[27] = fmod(x0, x1 + 1)
        # aux[28] = fract(x0) # Compile error
        aux[28] = heaviside(t - 0.5)
        aux[29] = gamma(x0 + 1)
        aux[30] = hypot(x0, x1)
        aux[31] = float(ilogb(x0 + 1))
        aux[32] = ldexp(x0, int(x1))
        aux[33] = lgamma(x0 + 1)
        aux[34] = log(x0 + 1)
        aux[35] = log1p(x0 + 1)
        aux[36] = log2(x0 + 1)
        aux[37] = log10(x0 + 1)
        aux[38] = 0.0  # logb(x0 + 1)
        aux[39] = 0.0  # mad(x0, x1, x2)
        aux[40] = 0.0  # nan()
        aux[41] = nextafter(x0, x1)
        aux[42] = pow(x0, x1)
        aux[43] = pown(x0, p0)
        aux[44] = powr(x0 + 0.5, x1 + 0.2)
        aux[45] = remainder(x0, x1 + 1)
        aux[46] = rint(x0)
        aux[47] = rootn(x0, p1)
        aux[48] = rsqrt(x0 + 1)
        aux[49] = sin(x0)
        aux[50] = sinh(x0)
        aux[51] = sinpi(x0)
        aux[52] = sqrt(x0)
        aux[53] = tan(x0)
        aux[54] = tanh(x0)
        aux[55] = tanpi(x0)
        aux[56] = trunc(x0)

        dx[0] = 0.0
        dx[1] = 0.0

    aux = [f"aux{num}" for num in range(57)]

    variables = {"x": 0.3, "y": 1.5, "z": 2.5}
    parameters = {"p0": 2.8, "p1": 4.55}

    sim = TrajectorySimulator(
        variables=variables,
        parameters=parameters,
        rhs_equation=rhs,
        dt=0.5,
        aux=aux,
        t_span=(0, 1),
        stepper=clode.Stepper.euler,
        max_store=3,
        max_steps=3,
    )

    trajectory = sim.trajectory()

    vars_arr = list(variables.values())
    params_arr = list(parameters.values())
    for t_index, t in enumerate([0, 0.5, 1]):
        aux_arr = [0.0] * len(aux)
        dx = [0.0] * len(vars_arr)
        w: List[float] = []
        rhs(t, vars_arr, params_arr, dx, aux_arr, w)
        for i, auxi in enumerate(aux):
            cl_aux_arr = trajectory.aux[auxi]
            assert np.isclose(
                cl_aux_arr[t_index], aux_arr[i], atol=1e-7
            ), f"Assertion failed for {auxi} at time {t}, expected {aux_arr[i]}, got {cl_aux_arr[t_index]}"


if __name__ == "__main__":
    test_opencl_builtins()
