from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from clode.runtime import CLDeviceType, CLVendor

from .errors import DoublePrecisionNotSupportedError, OpenCLDependencyError

try:
    import pyopencl as cl
except (ImportError, OSError):
    cl = None


def _require_opencl_binding() -> Any:
    if cl is None:
        raise OpenCLDependencyError(
            "pyopencl is required for clODE's OpenCL runtime. Install pyopencl into the active environment or reinstall clode with its default dependencies."
        )
    return cl


def _matches_device_type(device: Any, device_type: CLDeviceType | None) -> bool:
    if device_type is None:
        return True
    device_type_value = int(device_type)
    if device_type_value in {
        int(CLDeviceType.DEVICE_TYPE_DEFAULT),
        int(CLDeviceType.DEVICE_TYPE_ALL),
    }:
        return True
    return bool(int(device.type) & device_type_value)


def _matches_vendor(device: Any, vendor: CLVendor | None) -> bool:
    if vendor is None or int(vendor) == int(CLVendor.VENDOR_ANY):
        return True
    vendor_text = str(device.vendor).lower()
    vendor_map = {
        int(CLVendor.VENDOR_NVIDIA): ("nvidia",),
        int(CLVendor.VENDOR_AMD): ("amd", "advanced micro devices"),
        int(CLVendor.VENDOR_INTEL): ("intel",),
    }
    return any(token in vendor_text for token in vendor_map.get(int(vendor), ()))


@dataclass(slots=True)
class OpenCLRuntime:
    context: object
    queue: object
    platform: object
    device: object
    platform_id: int
    device_id: int
    program_cache: object = field(init=False, repr=False)
    struct_cache: dict[str, tuple[object, str]] = field(init=False, repr=False)

    def __post_init__(self) -> None:
        from .program_cache import ProgramCache

        self.program_cache = ProgramCache()
        self.struct_cache = {}

    @classmethod
    def create(
        cls,
        device_type: CLDeviceType | None = None,
        vendor: CLVendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
    ) -> OpenCLRuntime:
        opencl_binding = _require_opencl_binding()
        platforms = opencl_binding.get_platforms()

        if platform_id is not None:
            if device_type is not None:
                raise ValueError("Cannot specify device_type when platform_id is specified")
            if vendor is not None:
                raise ValueError("Cannot specify vendor when platform_id is specified")
            if device_id is None:
                raise ValueError("Must specify device_id when platform_id is specified")
            try:
                platform = platforms[platform_id]
            except IndexError as exc:
                raise ValueError(f"Invalid platform_id: {platform_id}") from exc
            devices = platform.get_devices()
            try:
                device = devices[device_id]
            except IndexError as exc:
                raise ValueError(f"Invalid device_id {device_id} for platform {platform_id}") from exc
            return cls._from_selected_device(platform, device, platform_id, device_id)

        if device_id is not None:
            raise ValueError("Must specify platform_id when specifying device_id")

        for selected_platform_id, platform in enumerate(platforms):
            for selected_device_id, device in enumerate(platform.get_devices()):
                if not _matches_device_type(device, device_type):
                    continue
                if not _matches_vendor(device, vendor):
                    continue
                return cls._from_selected_device(
                    platform,
                    device,
                    selected_platform_id,
                    selected_device_id,
                )

        raise ValueError("No matching OpenCL device found")

    @classmethod
    def _from_selected_device(
        cls,
        platform: object,
        device: object,
        platform_id: int,
        device_id: int,
    ) -> OpenCLRuntime:
        opencl_binding = _require_opencl_binding()
        context = opencl_binding.Context(devices=[device])
        queue = opencl_binding.CommandQueue(context)
        return cls(context, queue, platform, device, platform_id, device_id)

    def get_double_support(self) -> bool:
        extensions = str(getattr(self.device, "extensions", "")).lower().split()
        return (
            "cl_khr_fp64" in extensions
            or "cl_amd_fp64" in extensions
            or int(getattr(self.device, "double_fp_config", 0)) != 0
        )

    def describe(self) -> str:
        platform_name = str(getattr(self.platform, "name", "unknown platform"))
        device_name = str(getattr(self.device, "name", "unknown device"))
        return (
            f"platform_id={self.platform_id} ({platform_name}), "
            f"device_id={self.device_id} ({device_name})"
        )

    def get_max_memory_alloc_size(self) -> int:
        return int(self.device.max_mem_alloc_size)

    def get_device_cl_version(self) -> str:
        return str(self.device.version)

    def require_double_precision(self) -> None:
        if not self.get_double_support():
            raise DoublePrecisionNotSupportedError(str(self.device.name))
