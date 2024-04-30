# This is the package public interface (expose objects for "from clode import X")
# - keep only main simulation items here?
# - keep CL-related things in runtime?
from clode.runtime import (
    CLDeviceType,
    CLVendor,
    DeviceInfo,
    PlatformInfo,
    LogLevel,
    get_log_level,
    initialize_runtime,
    print_opencl,
    query_opencl,
    set_log_level,
    set_log_pattern,
)
from clode.cpp.clode_cpp_wrapper import (
    SolverParams,
    ObserverParams,
    ProblemInfo,
)
from clode.solver import Simulator, Stepper
from clode.features import FeatureSimulator, Observer, ObserverOutput
from clode.trajectory import TrajectorySimulator, TrajectoryOutput
from clode.function_converter import (
    OpenCLConverter,
    OpenCLRhsEquation,
    convert_str_to_opencl,
)
from clode.xpp_parser import convert_xpp_file, format_opencl_rhs, read_ode_parameters

# Import everything from opencl_builtins
from clode.opencl_builtins import (
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

__version__ = "0.8.0"

__all__ = [
    "CLDeviceType",
    "CLVendor",
    "DeviceInfo",
    "PlatformInfo",
    "initialize_runtime",
    "print_opencl",
    "query_opencl",
    "LogLevel",
    "get_log_level",
    "set_log_level",
    "set_log_pattern",
    "ProblemInfo",
    "SolverParams",
    "Stepper",
    "Simulator",
    "FeatureSimulator",
    "Observer",
    "ObserverParams",
    "ObserverOutput",
    "TrajectorySimulator",
    "TrajectoryOutput",
    "OpenCLConverter",
    "OpenCLRhsEquation",
    "convert_str_to_opencl",
    "convert_xpp_file",
    "format_opencl_rhs",
    "read_ode_parameters",
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
    "hypot",
    "heaviside",
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
