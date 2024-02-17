"""
CLODE C++/Python interface
"""
from __future__ import annotations
import typing
__all__ = ['CLDeviceType', 'CLVendor', 'DEVICE_TYPE_ACCELERATOR', 'DEVICE_TYPE_ALL', 'DEVICE_TYPE_CPU', 'DEVICE_TYPE_CUSTOM', 'DEVICE_TYPE_DEFAULT', 'DEVICE_TYPE_GPU', 'DeviceInfo', 'FeatureSimulatorBase', 'LogLevel', 'LoggerSingleton', 'ObserverParams', 'OpenCLResource', 'PlatformInfo', 'ProblemInfo', 'SimulatorBase', 'SolverParams', 'TrajectorySimulatorBase', 'VENDOR_AMD', 'VENDOR_ANY', 'VENDOR_INTEL', 'VENDOR_NVIDIA', 'critical', 'debug', 'err', 'get_logger', 'info', 'off', 'print_opencl', 'query_opencl', 'trace', 'warn']
class CLDeviceType:
    """
    Members:
    
      DEVICE_TYPE_ALL
    
      DEVICE_TYPE_CPU
    
      DEVICE_TYPE_GPU
    
      DEVICE_TYPE_ACCELERATOR
    
      DEVICE_TYPE_DEFAULT
    
      DEVICE_TYPE_CUSTOM
    """
    DEVICE_TYPE_ACCELERATOR: typing.ClassVar[CLDeviceType]  # value = <CLDeviceType.DEVICE_TYPE_ACCELERATOR: 8>
    DEVICE_TYPE_ALL: typing.ClassVar[CLDeviceType]  # value = <CLDeviceType.DEVICE_TYPE_ALL: 4294967295>
    DEVICE_TYPE_CPU: typing.ClassVar[CLDeviceType]  # value = <CLDeviceType.DEVICE_TYPE_CPU: 2>
    DEVICE_TYPE_CUSTOM: typing.ClassVar[CLDeviceType]  # value = <CLDeviceType.DEVICE_TYPE_CUSTOM: 16>
    DEVICE_TYPE_DEFAULT: typing.ClassVar[CLDeviceType]  # value = <CLDeviceType.DEVICE_TYPE_DEFAULT: 1>
    DEVICE_TYPE_GPU: typing.ClassVar[CLDeviceType]  # value = <CLDeviceType.DEVICE_TYPE_GPU: 4>
    __members__: typing.ClassVar[dict[str, CLDeviceType]]  # value = {'DEVICE_TYPE_ALL': <CLDeviceType.DEVICE_TYPE_ALL: 4294967295>, 'DEVICE_TYPE_CPU': <CLDeviceType.DEVICE_TYPE_CPU: 2>, 'DEVICE_TYPE_GPU': <CLDeviceType.DEVICE_TYPE_GPU: 4>, 'DEVICE_TYPE_ACCELERATOR': <CLDeviceType.DEVICE_TYPE_ACCELERATOR: 8>, 'DEVICE_TYPE_DEFAULT': <CLDeviceType.DEVICE_TYPE_DEFAULT: 1>, 'DEVICE_TYPE_CUSTOM': <CLDeviceType.DEVICE_TYPE_CUSTOM: 16>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class CLVendor:
    """
    Members:
    
      VENDOR_ANY
    
      VENDOR_NVIDIA
    
      VENDOR_AMD
    
      VENDOR_INTEL
    """
    VENDOR_AMD: typing.ClassVar[CLVendor]  # value = <CLVendor.VENDOR_AMD: 2>
    VENDOR_ANY: typing.ClassVar[CLVendor]  # value = <CLVendor.VENDOR_ANY: 0>
    VENDOR_INTEL: typing.ClassVar[CLVendor]  # value = <CLVendor.VENDOR_INTEL: 3>
    VENDOR_NVIDIA: typing.ClassVar[CLVendor]  # value = <CLVendor.VENDOR_NVIDIA: 1>
    __members__: typing.ClassVar[dict[str, CLVendor]]  # value = {'VENDOR_ANY': <CLVendor.VENDOR_ANY: 0>, 'VENDOR_NVIDIA': <CLVendor.VENDOR_NVIDIA: 1>, 'VENDOR_AMD': <CLVendor.VENDOR_AMD: 2>, 'VENDOR_INTEL': <CLVendor.VENDOR_INTEL: 3>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class DeviceInfo:
    compute_units: int
    device_available: int
    device_memory_size: int
    device_type: int
    device_type_str: str
    double_support: bool
    extensions: str
    max_clock: int
    max_memory_alloc_size: int
    max_work_group_size: int
    name: str
    vendor: str
    version: str
    def __repr__(self) -> str:
        """
        Device info string representation
        """
class FeatureSimulatorBase(SimulatorBase):
    def __init__(self, arg0: ProblemInfo, arg1: str, arg2: str, arg3: bool, arg4: OpenCLResource, arg5: str) -> None:
        ...
    def __repr__(self) -> str:
        """
        CLODEfeatures string representation
        """
    def build_cl(self) -> None:
        ...
    @typing.overload
    def features(self, arg0: bool) -> None:
        ...
    @typing.overload
    def features(self) -> None:
        ...
    def get_available_observers(self) -> list[str]:
        ...
    def get_f(self) -> list[float]:
        ...
    def get_feature_names(self) -> list[str]:
        ...
    def get_n_features(self) -> int:
        ...
    def get_observer_name(self) -> str:
        ...
    @typing.overload
    def initialize(self, arg0: list[float], arg1: list[float], arg2: list[float], arg3: SolverParams) -> None:
        """
        Initialize FeatureSimulatorBase
        """
    @typing.overload
    def initialize(self, arg0: list[float], arg1: list[float], arg2: list[float], arg3: SolverParams, arg4: ObserverParams) -> None:
        """
        Initialize FeatureSimulatorBase
        """
    def set_observer_params(self, arg0: ObserverParams) -> None:
        ...
class LogLevel:
    """
    Members:
    
      trace
    
      debug
    
      info
    
      warn
    
      err
    
      critical
    
      off
    """
    __members__: typing.ClassVar[dict[str, LogLevel]]  # value = {'trace': <LogLevel.trace: 0>, 'debug': <LogLevel.debug: 1>, 'info': <LogLevel.info: 2>, 'warn': <LogLevel.warn: 3>, 'err': <LogLevel.err: 4>, 'critical': <LogLevel.critical: 5>, 'off': <LogLevel.off: 6>}
    critical: typing.ClassVar[LogLevel]  # value = <LogLevel.critical: 5>
    debug: typing.ClassVar[LogLevel]  # value = <LogLevel.debug: 1>
    err: typing.ClassVar[LogLevel]  # value = <LogLevel.err: 4>
    info: typing.ClassVar[LogLevel]  # value = <LogLevel.info: 2>
    off: typing.ClassVar[LogLevel]  # value = <LogLevel.off: 6>
    trace: typing.ClassVar[LogLevel]  # value = <LogLevel.trace: 0>
    warn: typing.ClassVar[LogLevel]  # value = <LogLevel.warn: 3>
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class LoggerSingleton:
    def get_log_level(self) -> LogLevel:
        ...
    def set_log_level(self, arg0: LogLevel) -> None:
        ...
    def set_log_pattern(self, arg0: str) -> None:
        ...
class ObserverParams:
    e_var_ix: int
    f_var_ix: int
    maxEventCount: int
    def __init__(self, arg0: int, arg1: int, arg2: int, arg3: float, arg4: float, arg5: float, arg6: float, arg7: float, arg8: float, arg9: float, arg10: float) -> None:
        ...
    def __repr__(self) -> str:
        ...
class OpenCLResource:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: int) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: CLVendor) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: int, arg1: CLVendor) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: CLDeviceType, arg1: CLVendor) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: int, arg1: int) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: int, arg1: list[int]) -> None:
        ...
    def get_device_cl_version(self, arg0: int) -> str:
        """
        Get device CL version
        """
    def get_double_support(self, arg0: int) -> bool:
        """
        Get double support
        """
    def get_max_memory_alloc_size(self, arg0: int) -> int:
        """
        Get max memory alloc size
        """
    def print_devices(self) -> None:
        """
        Print device info to log
        """
class PlatformInfo:
    device_count: int
    device_info: list[DeviceInfo]
    name: str
    vendor: str
    version: str
    def __repr__(self) -> str:
        """
        Platform info string representation
        """
class ProblemInfo:
    aux: list[str]
    num_noise: int
    pars: list[str]
    src_file: str
    vars: list[str]
    @typing.overload
    def __init__(self, src_file: str, vars: list[str], pars: list[str], aux: list[str] = [], num_noise: int = 1) -> None:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        """
        clode Problem info string representation
        """
class SimulatorBase:
    def __init__(self, problem_info: ProblemInfo, stepper: str, cl_single_precision: bool, opencl_resource: OpenCLResource, clode_root: str) -> None:
        ...
    def build_cl(self) -> None:
        ...
    def get_available_steppers(self) -> list[str]:
        ...
    def get_problem_info(self) -> ProblemInfo:
        ...
    def get_program_string(self) -> str:
        ...
    def get_tspan(self) -> list[float]:
        ...
    def get_x0(self) -> list[float]:
        ...
    def get_xf(self) -> list[float]:
        ...
    def initialize(self, tspan: list[float], x0: list[float], pars: list[float], solver_params: SolverParams) -> None:
        """
        Initialize CLODE
        """
    def print_status(self) -> None:
        ...
    @typing.overload
    def seed_rng(self, seed: int) -> None:
        """
        Seed RNG
        """
    @typing.overload
    def seed_rng(self) -> None:
        """
        Seed RNG
        """
    def set_pars(self, arg0: list[float]) -> None:
        ...
    def set_problem_data(self, arg0: list[float], arg1: list[float]) -> None:
        ...
    def set_solver_params(self, arg0: SolverParams) -> None:
        ...
    def set_tspan(self, arg0: list[float]) -> None:
        ...
    def set_x0(self, arg0: list[float]) -> None:
        ...
    def shift_tspan(self) -> None:
        ...
    def shift_x0(self) -> None:
        ...
    def transient(self) -> None:
        ...
class SolverParams:
    abstol: float
    dt: float
    dtmax: float
    max_steps: int
    max_store: int
    nout: int
    reltol: float
    def __init__(self, dt: float = 0.1, dtmax: float = 0.5, abstol: float = 1e-06, reltol: float = 0.001, max_steps: int = 1000000, max_store: int = 1000000, nout: int = 1) -> None:
        ...
    def __repr__(self) -> str:
        """
        Solver params string representation
        """
class TrajectorySimulatorBase(SimulatorBase):
    def __init__(self, arg0: ProblemInfo, arg1: str, arg2: bool, arg3: OpenCLResource, arg4: str) -> None:
        ...
    def build_cl(self) -> None:
        ...
    def get_aux(self) -> list[float]:
        ...
    def get_dx(self) -> list[float]:
        ...
    def get_n_stored(self) -> list[int]:
        ...
    def get_t(self) -> list[float]:
        ...
    def get_x(self) -> list[float]:
        ...
    def initialize(self, arg0: list[float], arg1: list[float], arg2: list[float], arg3: SolverParams) -> None:
        """
        Initialize TrajectorySimulatorBase
        """
    def trajectory(self) -> None:
        ...
def get_logger() -> LoggerSingleton:
    """
    Get logger singleton instance
    """
def print_opencl() -> None:
    """
    Print OpenCL devices
    """
def query_opencl() -> list[PlatformInfo]:
    """
    Query OpenCL devices
    """
DEVICE_TYPE_ACCELERATOR: CLDeviceType  # value = <CLDeviceType.DEVICE_TYPE_ACCELERATOR: 8>
DEVICE_TYPE_ALL: CLDeviceType  # value = <CLDeviceType.DEVICE_TYPE_ALL: 4294967295>
DEVICE_TYPE_CPU: CLDeviceType  # value = <CLDeviceType.DEVICE_TYPE_CPU: 2>
DEVICE_TYPE_CUSTOM: CLDeviceType  # value = <CLDeviceType.DEVICE_TYPE_CUSTOM: 16>
DEVICE_TYPE_DEFAULT: CLDeviceType  # value = <CLDeviceType.DEVICE_TYPE_DEFAULT: 1>
DEVICE_TYPE_GPU: CLDeviceType  # value = <CLDeviceType.DEVICE_TYPE_GPU: 4>
VENDOR_AMD: CLVendor  # value = <CLVendor.VENDOR_AMD: 2>
VENDOR_ANY: CLVendor  # value = <CLVendor.VENDOR_ANY: 0>
VENDOR_INTEL: CLVendor  # value = <CLVendor.VENDOR_INTEL: 3>
VENDOR_NVIDIA: CLVendor  # value = <CLVendor.VENDOR_NVIDIA: 1>
critical: LogLevel  # value = <LogLevel.critical: 5>
debug: LogLevel  # value = <LogLevel.debug: 1>
err: LogLevel  # value = <LogLevel.err: 4>
info: LogLevel  # value = <LogLevel.info: 2>
off: LogLevel  # value = <LogLevel.off: 6>
trace: LogLevel  # value = <LogLevel.trace: 0>
warn: LogLevel  # value = <LogLevel.warn: 3>
