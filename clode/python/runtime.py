import os

from . import clode_cpp_wrapper as _clode

_runtime = None
_clode_root_dir: str = os.path.join(os.path.dirname(__file__), "cpp", "")


def get_runtime(
    device_type: _clode.cl_device_type = _clode.cl_device_type.DEVICE_TYPE_DEFAULT,
    vendor: _clode.cl_vendor = _clode.cl_vendor.VENDOR_ANY,
) -> _clode.opencl_resource:
    return initialise_runtime(device_type, vendor)


def initialise_runtime(
    device_type: _clode.cl_device_type = _clode.cl_device_type.DEVICE_TYPE_DEFAULT,
    vendor: _clode.cl_vendor = _clode.cl_vendor.VENDOR_ANY,
) -> _clode.opencl_resource:
    global _runtime
    if _runtime is None:
        _runtime = _clode.opencl_resource(device_type, vendor)
    return _runtime


def initialise_runtime_by_platform_id(
    platformID: int, deviceID: int
) -> _clode.opencl_resource:
    global _runtime
    if _runtime is None:
        _runtime = _clode.opencl_resource(platformID, deviceID)
    return _runtime


def initialise_runtime_by_device_ids(
    platformID: int, device_ids: list[int]
) -> _clode.opencl_resource:
    global _runtime
    if _runtime is None:
        _runtime = _clode.opencl_resource(platformID, device_ids)
    return _runtime


def reset_runtime():
    global _runtime
    _runtime = None


def get_cpp():
    return _clode


def print_devices():
    return get_runtime().print()
