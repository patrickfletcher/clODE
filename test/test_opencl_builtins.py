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
    fma,
    fmod,
    fract,
    gamma,
    hypot,
    ilogb,
    ldexp,
    lgamma,
    log,
    log1p,
    log2,
    log10,
    logb,
    mad,
    nan,
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


def test_opencl_builtins():
    def rhs(
        t: float,
        x: list[float],
        p: list[float],
        dx: list[float],
        aux: list[float],
        w: list[float],
    ) -> None:
        #  Test every OpenCL builtin function and place the result in aux

        p0: int = int(p[0])
        p1: int = int(p[1])
        aux[0] = acos(x[0])
        aux[1] = acosh(x[0])
        aux[2] = acospi(x[0])
        aux[3] = asin(x[0])
        aux[4] = asinh(x[0])
        aux[5] = asinpi(x[0])
        aux[6] = atan(x[0])
        aux[7] = atan2(x[0], x[1])
        aux[8] = atan2pi(x[0], x[1])
        aux[9] = atanh(x[0])
        aux[10] = atanpi(x[0])
        aux[11] = cbrt(x[0])
        aux[12] = ceil(x[0])
        aux[13] = copysign(x[0], x[1])
        aux[14] = cos(x[0])
        aux[15] = cosh(x[0])
        aux[16] = cospi(x[0])
        aux[17] = erf(x[0])
        aux[18] = erfc(x[0])
        aux[19] = exp(x[0])
        aux[20] = exp2(x[0])
        aux[21] = exp10(x[0])
        aux[22] = expm1(x[0])
        aux[23] = fabs(x[0])
        aux[24] = fdim(x[0], x[1])
        aux[25] = floor(x[0])
        aux[26] = fma(x[0], x[1], x[2])
        aux[27] = fmod(x[0], x[1])
        # aux[28] = fract(x[0])
        aux[29] = gamma(x[0])
        aux[30] = hypot(x[0], x[1])
        aux[31] = float(ilogb(x[0]))
        aux[32] = ldexp(x[0], x[1])
        aux[33] = lgamma(x[0])
        aux[34] = log(x[0])
        aux[35] = log1p(x[0])
        aux[36] = log2(x[0])
        aux[37] = log10(x[0])
        aux[38] = logb(x[0])
        aux[39] = mad(x[0], x[1], x[2])
        # aux[40] = nan()
        aux[41] = nextafter(x[0], x[1])
        aux[42] = pow(x[0], x[1])
        aux[43] = pown(x[0], p0)
        aux[44] = powr(x[0], x[1])
        aux[45] = remainder(x[0], x[1])
        aux[46] = rint(x[0])
        aux[47] = rootn(x[0], p1)
        aux[48] = rsqrt(x[0])
        aux[49] = sin(x[0])
        aux[50] = sinh(x[0])
        aux[51] = sinpi(x[0])
        aux[52] = sqrt(x[0])
        aux[53] = tan(x[0])
        aux[54] = tanh(x[0])
        aux[55] = tanpi(x[0])
        aux[56] = trunc(x[0])

        dx[0] = 0.0
        dx[1] = 0.0

    aux = [f"aux{num}" for num in range(57)]

    sim = TrajectorySimulator(
        variables={"x": 0.3, "y": 1.5},
        parameters={"p0": 2.8, "p1": 4.55},
        rhs_equation=rhs,
        dt=0.5,
        aux=aux,
        t_span=(0, 1),
    )

    sim.trajectory()
    state = sim.get_final_state()
    aux_values = sim.get_aux()
    print(0)
