from __future__ import annotations

import os

from clode.cpp import clode_cpp_wrapper as _clode
from clode.cpp.clode_cpp_wrapper import (
    CLDeviceType,
    CLVendor,
    DeviceInfo,
    PlatformInfo,
    LogLevel,
    ProblemInfo,
    SolverParams,
    ObserverParams,
    OpenCLResource,
    SimulatorBase,
    FeatureSimulatorBase,
    TrajectorySimulatorBase,
    print_opencl,
    query_opencl,
    get_logger,
)

_clode_root_dir: str = os.path.join(os.path.dirname(__file__), "cpp", "")

DEFAULT_LOG_LEVEL = LogLevel.warn

def initialize_runtime(
    device_type: CLDeviceType | None,
    vendor: CLVendor | None,
    platform_id: int | None,
    device_id: int | None,
    device_ids: list[int] | None,
) -> OpenCLResource:
    if platform_id is not None:
        if device_type is not None:
            raise ValueError("Cannot specify device_type when platform_id is specified")
        if vendor is not None:
            raise ValueError("Cannot specify vendor when platform_id is specified")
        if device_id is not None and device_ids is not None:
            raise ValueError("Cannot specify both device_id and device_ids")
        if device_id is None and device_ids is None:
            raise ValueError("Must specify one of device_id and device_ids")
        if device_id is not None:
            return OpenCLResource(platform_id, device_id)
        if device_ids is not None:
            return OpenCLResource(platform_id, device_ids)
        raise ValueError("Must specify one of device_id and device_ids")
    elif device_id is not None:
        raise ValueError("Must specify platform_id when specifying device_id")
    elif device_ids is not None:
        raise ValueError("Must specify platform_id when specifying device_ids")
    else:
        if device_type is None:
            device_type = CLDeviceType.DEVICE_TYPE_DEFAULT
        if vendor is None:
            vendor = CLVendor.VENDOR_ANY
        return OpenCLResource(device_type, vendor)




def get_log_level() -> LogLevel:
    return get_logger().get_log_level()


def set_log_level(level: LogLevel) -> None:
    get_logger().set_log_level(level)


def set_log_pattern(pattern: str) -> None:
    get_logger().set_log_pattern(pattern)


set_log_level(DEFAULT_LOG_LEVEL)

__all__ = [
    "CLDeviceType",
    "CLVendor",
    "DeviceInfo",
    "PlatformInfo",
    "DEFAULT_LOG_LEVEL",
    "get_log_level",
    "initialize_runtime",
    "LogLevel",
    "set_log_level",
    "set_log_pattern",
    "ProblemInfo",
    "SolverParams",
    "ObserverParams",
    "SimulatorBase",
    "OpenCLResource",
    "FeatureSimulatorBase",
    "TrajectorySimulatorBase",
    "print_opencl",
    "query_opencl",
]
