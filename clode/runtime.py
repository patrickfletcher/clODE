from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import os
from typing import Any, Sequence

_clode_root_dir: str = os.path.join(os.path.dirname(__file__), "kernels", "")
_BACKEND_ENVVAR = "_CLODE_BACKEND"
_DEFAULT_BACKEND = "pyopencl"


class CLDeviceType(IntEnum):
    DEVICE_TYPE_ALL = -1
    DEVICE_TYPE_DEFAULT = 1
    DEVICE_TYPE_CPU = 2
    DEVICE_TYPE_GPU = 4
    DEVICE_TYPE_ACCELERATOR = 8
    DEVICE_TYPE_CUSTOM = 16


class CLVendor(IntEnum):
    VENDOR_ANY = 0
    VENDOR_NVIDIA = 1
    VENDOR_AMD = 2
    VENDOR_INTEL = 3


class LogLevel(IntEnum):
    trace = 0
    debug = 1
    info = 2
    warn = 3
    err = 4
    critical = 5
    off = 6


@dataclass(slots=True)
class DeviceInfo:
    name: str
    vendor: str
    version: str
    device_type: int
    device_type_str: str
    compute_units: int
    max_clock: int
    max_work_group_size: int
    device_memory_size: int
    max_memory_alloc_size: int
    extensions: str
    double_support: bool
    device_available: bool

    def __repr__(self) -> str:
        return (
            "<device_info("
            f"name={self.name}, "
            f"vendor={self.vendor}, "
            f"version={self.version}, "
            f"device_type={self.device_type_str}, "
            f"compute_units={self.compute_units}, "
            f"max_clock={self.max_clock}, "
            f"max_work_group_size={self.max_work_group_size}, "
            f"device_memory_size={self.device_memory_size}, "
            f"max_memory_alloc_size={self.max_memory_alloc_size}, "
            f"extensions={self.extensions}, "
            f"double_support={int(self.double_support)}, "
            f"device_available={int(self.device_available)}"
            ")>"
        )


@dataclass(slots=True)
class PlatformInfo:
    name: str
    vendor: str
    version: str
    device_count: int
    device_info: list[DeviceInfo]

    def __repr__(self) -> str:
        return (
            "<platform_info("
            f"name={self.name}, "
            f"vendor={self.vendor}, "
            f"version={self.version}, "
            f"device_count={self.device_count}"
            ")>"
        )

DEFAULT_LOG_LEVEL = LogLevel.warn
_LOG_LEVEL = DEFAULT_LOG_LEVEL
_LOG_PATTERN: str | None = None

_DEVICE_TYPE_VALUES = {member.value for member in CLDeviceType}
_VENDOR_VALUES = {member.value for member in CLVendor}


def resolve_backend_name(requested_backend: str | None = None) -> str:
    backend_name = (
        requested_backend
        if requested_backend is not None
        else os.getenv(_BACKEND_ENVVAR, _DEFAULT_BACKEND)
    )
    if backend_name != "pyopencl":
        raise ValueError(
            f"Unsupported clODE backend '{backend_name}'. Supported backends: ['pyopencl']"
        )
    return backend_name


def _load_pyopencl() -> Any | None:
    try:
        import pyopencl as cl
    except (ImportError, OSError):
        return None
    return cl


def _require_pyopencl() -> Any:
    pyopencl = _load_pyopencl()
    if pyopencl is None:
        from ._pyopencl.errors import PyOpenCLDependencyError

        raise PyOpenCLDependencyError(
            "pyopencl is required for this clODE installation. Install pyopencl into the active environment or reinstall the package with its default dependencies."
        )
    return pyopencl


def _coerce_device_type(device_type: CLDeviceType | int | None) -> CLDeviceType | None:
    if device_type is None:
        return None
    return CLDeviceType(int(device_type))


def _coerce_vendor(vendor: CLVendor | int | None) -> CLVendor | None:
    if vendor is None:
        return None
    return CLVendor(int(vendor))


def _coerce_log_level(level: LogLevel | int) -> LogLevel:
    return LogLevel(int(level))


def _normalize_runtime_selection(
    device_type: CLDeviceType | int | None,
    vendor: CLVendor | int | None,
    platform_id: int | None,
    device_id: int | None,
    device_ids: Sequence[int] | None,
) -> tuple[CLDeviceType | None, CLVendor | None, int | None, int | None, tuple[int, ...] | None]:
    normalized_device_type = _coerce_device_type(device_type)
    normalized_vendor = _coerce_vendor(vendor)
    normalized_platform_id = None if platform_id is None else int(platform_id)
    normalized_device_id = None if device_id is None else int(device_id)
    normalized_device_ids = (
        None if device_ids is None else tuple(int(candidate) for candidate in device_ids)
    )

    if normalized_platform_id is not None:
        if normalized_device_type is not None:
            raise ValueError("Cannot specify device_type when platform_id is specified")
        if normalized_vendor is not None:
            raise ValueError("Cannot specify vendor when platform_id is specified")
        if normalized_device_id is not None and normalized_device_ids is not None:
            raise ValueError("Cannot specify both device_id and device_ids")
        if normalized_device_id is None and normalized_device_ids is None:
            raise ValueError("Must specify one of device_id and device_ids")
    elif normalized_device_id is not None:
        raise ValueError("Must specify platform_id when specifying device_id")
    elif normalized_device_ids is not None:
        raise ValueError("Must specify platform_id when specifying device_ids")

    return (
        normalized_device_type,
        normalized_vendor,
        normalized_platform_id,
        normalized_device_id,
        normalized_device_ids,
    )


def _parse_opencl_resource_args(
    args: tuple[Any, ...]
) -> tuple[CLDeviceType | None, CLVendor | None, int | None, int | None, tuple[int, ...] | None]:
    if len(args) == 0:
        return None, None, None, None, None

    if len(args) == 1:
        arg0 = args[0]
        if isinstance(arg0, CLDeviceType):
            return arg0, None, None, None, None
        if isinstance(arg0, CLVendor):
            return None, arg0, None, None, None
        if isinstance(arg0, int):
            matches_device_type = int(arg0) in _DEVICE_TYPE_VALUES
            matches_vendor = int(arg0) in _VENDOR_VALUES
            if matches_device_type and not matches_vendor:
                return CLDeviceType(int(arg0)), None, None, None, None
            if matches_vendor and not matches_device_type:
                return None, CLVendor(int(arg0)), None, None, None
            raise TypeError(
                "Ambiguous single integer argument for OpenCLResource; use CLDeviceType or CLVendor explicitly."
            )
        raise TypeError("Unsupported OpenCLResource constructor arguments")

    if len(args) == 2:
        arg0, arg1 = args
        if isinstance(arg1, Sequence) and not isinstance(arg1, (str, bytes)):
            return None, None, int(arg0), None, tuple(int(candidate) for candidate in arg1)
        if isinstance(arg1, CLVendor):
            return _coerce_device_type(arg0), arg1, None, None, None
        if isinstance(arg0, int) and isinstance(arg1, int):
            return None, None, int(arg0), int(arg1), None
        raise TypeError("Unsupported OpenCLResource constructor arguments")

    raise TypeError("OpenCLResource accepts at most two positional arguments")


def _double_support_from_device(device: Any) -> bool:
    extensions = str(getattr(device, "extensions", "")).lower().split()
    return (
        "cl_khr_fp64" in extensions
        or "cl_amd_fp64" in extensions
        or int(getattr(device, "double_fp_config", 0)) != 0
    )


def _device_type_str(device_type: int) -> str:
    if device_type == int(CLDeviceType.DEVICE_TYPE_ALL):
        return "ALL"
    parts: list[str] = []
    if device_type & int(CLDeviceType.DEVICE_TYPE_DEFAULT):
        parts.append("DEFAULT")
    if device_type & int(CLDeviceType.DEVICE_TYPE_CPU):
        parts.append("CPU")
    if device_type & int(CLDeviceType.DEVICE_TYPE_GPU):
        parts.append("GPU")
    if device_type & int(CLDeviceType.DEVICE_TYPE_ACCELERATOR):
        parts.append("ACCELERATOR")
    if device_type & int(CLDeviceType.DEVICE_TYPE_CUSTOM):
        parts.append("CUSTOM")
    return "|".join(parts) if parts else str(device_type)


def _device_info_from_cpp(device_info: Any) -> DeviceInfo:
    return DeviceInfo(
        name=str(device_info.name),
        vendor=str(device_info.vendor),
        version=str(device_info.version),
        device_type=int(device_info.device_type),
        device_type_str=str(device_info.device_type_str),
        compute_units=int(device_info.compute_units),
        max_clock=int(device_info.max_clock),
        max_work_group_size=int(device_info.max_work_group_size),
        device_memory_size=int(device_info.device_memory_size),
        max_memory_alloc_size=int(device_info.max_memory_alloc_size),
        extensions=str(device_info.extensions),
        double_support=bool(device_info.double_support),
        device_available=bool(device_info.device_available),
    )


def _query_opencl_pyopencl() -> list[PlatformInfo]:
    pyopencl = _require_pyopencl()
    platforms: list[PlatformInfo] = []
    for platform in pyopencl.get_platforms():
        device_info: list[DeviceInfo] = []
        for device in platform.get_devices():
            device_type = int(getattr(device, "type", 0))
            device_info.append(
                DeviceInfo(
                    name=str(getattr(device, "name", "unknown device")),
                    vendor=str(getattr(device, "vendor", "unknown vendor")),
                    version=str(getattr(device, "version", "unknown version")),
                    device_type=device_type,
                    device_type_str=_device_type_str(device_type),
                    compute_units=int(getattr(device, "max_compute_units", 0)),
                    max_clock=int(getattr(device, "max_clock_frequency", 0)),
                    max_work_group_size=int(getattr(device, "max_work_group_size", 0)),
                    device_memory_size=int(getattr(device, "global_mem_size", 0)),
                    max_memory_alloc_size=int(getattr(device, "max_mem_alloc_size", 0)),
                    extensions=str(getattr(device, "extensions", "")),
                    double_support=_double_support_from_device(device),
                    device_available=bool(getattr(device, "available", False)),
                )
            )
        platforms.append(
            PlatformInfo(
                name=str(getattr(platform, "name", "unknown platform")),
                vendor=str(getattr(platform, "vendor", "unknown vendor")),
                version=str(getattr(platform, "version", "unknown version")),
                device_count=len(device_info),
                device_info=device_info,
            )
        )
    return platforms


def _emit_opencl_report(platforms: list[PlatformInfo]) -> None:
    print("Querying OpenCL platforms...")
    print(f"Number of platforms found: {len(platforms)}")
    for platform_index, platform in enumerate(platforms):
        print()
        print(f"Platform {platform_index}. ------------------------------")
        print(f"Name:    {platform.name}")
        print(f"Vendor:  {platform.vendor}")
        print(f"Version: {platform.version}")
        for device_index, device in enumerate(platform.device_info):
            memory_mb = device.device_memory_size // (1024 * 1024)
            max_alloc_mb = device.max_memory_alloc_size // (1024 * 1024)
            print()
            print(f"Device {device_index}. --------------------")
            print(f"Name:   {device.name}")
            print(f"Type:   {device.device_type_str}")
            print(f"Vendor: {device.vendor}")
            print(f"Version: {device.version}")
            print(f"Compute units (CUs): {device.compute_units}")
            print(f"Clock frequency:     {device.max_clock} MHz")
            print(f"Global memory size:  {memory_mb} MB")
            print(f"Max allocation size: {max_alloc_mb} MB")
            print(f"Max work group/CU:   {device.max_work_group_size}")
            print(f"Double support:      {str(device.double_support).lower()}")
            print(f"Device available:    {str(device.device_available).lower()}")


class OpenCLResource:
    def __init__(self, *args: Any) -> None:
        (
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        ) = _parse_opencl_resource_args(args)
        self._initialize(device_type, vendor, platform_id, device_id, device_ids)

    @classmethod
    def from_selection(
        cls,
        device_type: CLDeviceType | int | None,
        vendor: CLVendor | int | None,
        platform_id: int | None,
        device_id: int | None,
        device_ids: Sequence[int] | None,
    ) -> OpenCLResource:
        resource = cls.__new__(cls)
        resource._initialize(device_type, vendor, platform_id, device_id, device_ids)
        return resource

    def _initialize(
        self,
        device_type: CLDeviceType | int | None,
        vendor: CLVendor | int | None,
        platform_id: int | None,
        device_id: int | None,
        device_ids: Sequence[int] | None,
    ) -> None:
        (
            self._device_type,
            self._vendor,
            self._platform_id,
            self._device_id,
            self._device_ids,
        ) = _normalize_runtime_selection(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )
        resolve_backend_name()
        self._pyopencl_runtime: Any | None = None
        self._initialize_pyopencl_runtime()

    def _initialize_pyopencl_runtime(self) -> None:
        from ._pyopencl.runtime import OpenCLRuntime

        if self._device_ids is not None:
            if len(self._device_ids) != 1:
                raise ValueError(
                    "PyOpenCL backend currently supports exactly one selected device"
                )
            self._device_id = self._device_ids[0]

        if self._platform_id is not None:
            self._pyopencl_runtime = OpenCLRuntime.create(
                platform_id=self._platform_id,
                device_id=self._device_id,
            )
            return

        self._pyopencl_runtime = OpenCLRuntime.create(
            device_type=(
                CLDeviceType.DEVICE_TYPE_DEFAULT
                if self._device_type is None
                else self._device_type
            ),
            vendor=CLVendor.VENDOR_ANY if self._vendor is None else self._vendor,
        )
        self._platform_id = self._pyopencl_runtime.platform_id
        self._device_id = self._pyopencl_runtime.device_id

    def _ensure_selected_pyopencl_device(self, device_id: int) -> None:
        if self._pyopencl_runtime is None:
            return
        selected_device_id = 0 if self._device_id is None else self._device_id
        if device_id not in {0, selected_device_id}:
            raise ValueError(
                "PyOpenCL runtime resource is bound to a single selected device"
            )

    def get_device_cl_version(self, device_id: int) -> str:
        self._ensure_selected_pyopencl_device(device_id)
        return self._pyopencl_runtime.get_device_cl_version()

    def get_double_support(self, device_id: int) -> bool:
        self._ensure_selected_pyopencl_device(device_id)
        return self._pyopencl_runtime.get_double_support()

    def get_max_memory_alloc_size(self, device_id: int) -> int:
        self._ensure_selected_pyopencl_device(device_id)
        return self._pyopencl_runtime.get_max_memory_alloc_size()

    def print_devices(self) -> None:
        if get_log_level() == LogLevel.off:
            return
        _emit_opencl_report(_query_opencl_pyopencl())


def initialize_runtime(
    device_type: CLDeviceType | None,
    vendor: CLVendor | None,
    platform_id: int | None,
    device_id: int | None,
    device_ids: list[int] | None,
) -> OpenCLResource:
    return OpenCLResource.from_selection(
        device_type,
        vendor,
        platform_id,
        device_id,
        device_ids,
    )


def get_log_level() -> LogLevel:
    return _LOG_LEVEL


def set_log_level(level: LogLevel) -> None:
    global _LOG_LEVEL

    _LOG_LEVEL = _coerce_log_level(level)


def set_log_pattern(pattern: str) -> None:
    global _LOG_PATTERN

    _LOG_PATTERN = pattern


def query_opencl() -> list[PlatformInfo]:
    resolve_backend_name()
    return _query_opencl_pyopencl()


def print_opencl() -> None:
    old_level = get_log_level()
    if old_level == LogLevel.off:
        return
    if old_level > LogLevel.info:
        set_log_level(LogLevel.info)
    try:
        resolve_backend_name()
        _emit_opencl_report(_query_opencl_pyopencl())
    finally:
        if old_level > LogLevel.info:
            set_log_level(old_level)


__all__ = [
    "CLDeviceType",
    "CLVendor",
    "DeviceInfo",
    "PlatformInfo",
    "OpenCLResource",
    "initialize_runtime",
    "print_opencl",
    "query_opencl",
    "DEFAULT_LOG_LEVEL",
    "LogLevel",
    "set_log_level",
    "set_log_pattern",
    "get_log_level",
]
