import os

from . import clode_cpp_wrapper as _clode  # type: ignore
from .clode_cpp_wrapper import log_level, print_opencl, query_opencl  # type: ignore

_clode_root_dir: str = os.path.join(os.path.dirname(__file__), "cpp", "")

DEFAULT_LOG_LEVEL = _clode.log_level.warn


def initialise_runtime(
    device_type: _clode.cl_device_type | None,
    vendor: _clode.cl_vendor | None,
    platform_id: int | None,
    device_id: int | None,
    device_ids: list[int] | None,
) -> _clode.opencl_resource:
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
            return _clode.opencl_resource(platform_id, device_id)
        if device_ids is not None:
            return _clode.opencl_resource(platform_id, device_ids)
    elif device_id is not None:
        raise ValueError("Must specify platform_id when specifying device_id")
    elif device_ids is not None:
        raise ValueError("Must specify platform_id when specifying device_ids")
    else:
        if device_type is None:
            device_type = _clode.cl_device_type.DEVICE_TYPE_DEFAULT
        if vendor is None:
            vendor = _clode.cl_vendor.VENDOR_ANY
        return _clode.opencl_resource(device_type, vendor)


def get_cpp():
    return _clode


cl_device_type = _clode.cl_device_type
cl_vendor = _clode.cl_vendor

log_level = _clode.log_level


def get_log_level() -> log_level:
    return _clode.get_logger().get_log_level()


def set_log_level(level: log_level):
    return _clode.get_logger().set_log_level(level)


def set_log_pattern(pattern: str):
    return _clode.get_logger().set_log_pattern(pattern)


set_log_level(DEFAULT_LOG_LEVEL)
