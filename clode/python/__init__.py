from . import clode_cpp_wrapper as _clode  # type: ignore
from .clode_cpp_wrapper import print_opencl, query_opencl  # type: ignore
from .features import CLODEFeatures
from .observer import Observer, ObserverOutput
from .problem_info import ProblemInfo
from .runtime import (
    get_cpp,
    get_runtime,
    initialise_runtime,
    initialise_runtime_by_device_ids,
    initialise_runtime_by_platform_id,
    print_devices,
    reset_runtime,
    cl_device_type,
    cl_vendor,
)
from .stepper import Stepper
from .trajectory import CLODETrajectory
from .xpp_parser import convert_xpp_file, format_opencl_rhs, read_ode_parameters

log_level = _clode.log_level
set_log_level = _clode.set_log_level
set_log_pattern = _clode.set_log_pattern

__version__ = "0.4.0"
