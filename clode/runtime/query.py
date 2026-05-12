from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .logging import LogLevel, get_log_level, set_log_level


@dataclass(slots=True)
class DeviceInfo:
    """Summary of one OpenCL device visible to the active runtime."""

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
    """Summary of one OpenCL platform and its visible devices."""

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


def _load_opencl_binding() -> Any | None:
    try:
        import pyopencl as opencl_binding
    except (ImportError, OSError):
        return None
    return opencl_binding


def _require_opencl_binding() -> Any:
    opencl_binding = _load_opencl_binding()
    if opencl_binding is None:
        from .._opencl.errors import OpenCLDependencyError

        raise OpenCLDependencyError(
            "pyopencl is required for this clODE installation. Install pyopencl into the active environment or reinstall the package with its default dependencies."
        )
    return opencl_binding


def _double_support_from_device(device: Any) -> bool:
    extensions = str(getattr(device, "extensions", "")).lower().split()
    return (
        "cl_khr_fp64" in extensions
        or "cl_amd_fp64" in extensions
        or int(getattr(device, "double_fp_config", 0)) != 0
    )


def _device_type_str(device_type: int) -> str:
    from .selection import CLDeviceType

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


def _query_opencl_runtime() -> list[PlatformInfo]:
    opencl_binding = _require_opencl_binding()
    platforms: list[PlatformInfo] = []
    for platform in opencl_binding.get_platforms():
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


def query_opencl() -> list[PlatformInfo]:
    """Return the OpenCL platforms and devices visible to clODE."""

    return _query_opencl_runtime()


def print_opencl() -> None:
    """Print a human-readable report of visible OpenCL platforms and devices."""

    old_level = get_log_level()
    if old_level == LogLevel.off:
        return
    if old_level > LogLevel.info:
        set_log_level(LogLevel.info)
    try:
        _emit_opencl_report(_query_opencl_runtime())
    finally:
        if old_level > LogLevel.info:
            set_log_level(old_level)


__all__ = [
    "DeviceInfo",
    "PlatformInfo",
    "print_opencl",
    "query_opencl",
]