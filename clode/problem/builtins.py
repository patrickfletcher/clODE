from math import (
    acos,
    acosh,
    asin,
    asinh,
    atan,
    atan2,
    atanh,
    ceil,
    copysign,
    cos,
    cosh,
    erf,
    erfc,
    exp,
    expm1,
    fabs,
    floor,
    fmod,
    gamma,
    hypot,
    ldexp,
    lgamma,
    log,
    log1p,
    log2,
    log10,
    pi,
    pow,
    sin,
    sinh,
    sqrt,
    tan,
    tanh,
    trunc,
)

from numpy import cbrt, exp2, remainder


def acospi(x: float) -> float:
    return acos(x) / pi


def asinpi(x: float) -> float:
    return asin(x) / pi


def atanpi(x: float) -> float:
    return atan(x) / pi


def atan2pi(y: float, x: float) -> float:
    return atan2(y, x) / pi


def cospi(x: float) -> float:
    return cos(x * pi)


def exp10(x: float) -> float:
    return pow(10, x)


def fdim(x: float, y: float) -> float:
    return max(x - y, 0)


def ilogb(x: float) -> int:
    return int(log2(x))


def heaviside(x: float) -> float:
    return 1.0 if x >= 0.0 else 0.0


def nextafter(x: float, y: float) -> float:
    return x + copysign(0, y - x)


def pown(x: float, y: int) -> float:
    return pow(x, y)


def powr(x: float, y: float) -> float:
    return pow(x, y)


def rint(x: float) -> float:
    return round(x)


def rootn(x: float, y: int) -> float:
    return pow(x, 1 / y)


def rsqrt(x: float) -> float:
    return 1 / sqrt(x)


def sinpi(x: float) -> float:
    return sin(x * pi)


def tanpi(x: float) -> float:
    return tan(x * pi)


__all__ = [
    "acos",
    "acosh",
    "acospi",
    "asin",
    "asinh",
    "asinpi",
    "atan",
    "atan2",
    "atan2pi",
    "atanh",
    "atanpi",
    "cbrt",
    "ceil",
    "copysign",
    "cos",
    "cosh",
    "cospi",
    "erf",
    "erfc",
    "exp",
    "exp2",
    "exp10",
    "expm1",
    "fabs",
    "fdim",
    "floor",
    "fmod",
    "gamma",
    "heaviside",
    "hypot",
    "ilogb",
    "ldexp",
    "lgamma",
    "log",
    "log1p",
    "log2",
    "log10",
    "nextafter",
    "pow",
    "pown",
    "powr",
    "remainder",
    "rint",
    "rootn",
    "rsqrt",
    "sin",
    "sinh",
    "sinpi",
    "sqrt",
    "tan",
    "tanh",
    "tanpi",
    "trunc",
]