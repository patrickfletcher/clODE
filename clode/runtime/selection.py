from __future__ import annotations

from dataclasses import dataclass
from enum import IntEnum
import os
from typing import Any, Sequence

from .logging import LogLevel, get_log_level
from .query import _emit_opencl_report, _query_opencl_runtime

_clode_root_dir: str = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "kernels", ""
)


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


@dataclass(frozen=True, slots=True)
class RuntimeSelection:
    device_type: CLDeviceType | None
    vendor: CLVendor | None
    platform_id: int | None
    device_id: int | None
    device_ids: tuple[int, ...] | None


_DEVICE_TYPE_VALUES = {member.value for member in CLDeviceType}
_VENDOR_VALUES = {member.value for member in CLVendor}


def _coerce_device_type(device_type: CLDeviceType | int | None) -> CLDeviceType | None:
    if device_type is None:
        return None
    return CLDeviceType(int(device_type))


def _coerce_vendor(vendor: CLVendor | int | None) -> CLVendor | None:
    if vendor is None:
        return None
    return CLVendor(int(vendor))


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
        self._opencl_runtime: Any | None = None
        self._initialize_opencl_runtime()

    def _initialize_opencl_runtime(self) -> None:
        from .._opencl.runtime import OpenCLRuntime

        if self._device_ids is not None:
            if len(self._device_ids) != 1:
                raise ValueError(
                    "OpenCL execution currently supports exactly one selected device"
                )
            self._device_id = self._device_ids[0]

        if self._platform_id is not None:
            self._opencl_runtime = OpenCLRuntime.create(
                platform_id=self._platform_id,
                device_id=self._device_id,
            )
            return

        self._opencl_runtime = OpenCLRuntime.create(
            device_type=(
                CLDeviceType.DEVICE_TYPE_DEFAULT
                if self._device_type is None
                else self._device_type
            ),
            vendor=CLVendor.VENDOR_ANY if self._vendor is None else self._vendor,
        )
        self._platform_id = self._opencl_runtime.platform_id
        self._device_id = self._opencl_runtime.device_id

    def _ensure_selected_opencl_device(self, device_id: int) -> None:
        if self._opencl_runtime is None:
            return
        selected_device_id = 0 if self._device_id is None else self._device_id
        if device_id not in {0, selected_device_id}:
            raise ValueError(
                "OpenCL runtime resource is bound to a single selected device"
            )

    def get_device_cl_version(self, device_id: int) -> str:
        self._ensure_selected_opencl_device(device_id)
        return self._opencl_runtime.get_device_cl_version()

    def get_double_support(self, device_id: int) -> bool:
        self._ensure_selected_opencl_device(device_id)
        return self._opencl_runtime.get_double_support()

    def get_max_memory_alloc_size(self, device_id: int) -> int:
        self._ensure_selected_opencl_device(device_id)
        return self._opencl_runtime.get_max_memory_alloc_size()

    def print_devices(self) -> None:
        if get_log_level() == LogLevel.off:
            return
        _emit_opencl_report(_query_opencl_runtime())


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


__all__ = [
    "CLDeviceType",
    "CLVendor",
    "OpenCLResource",
    "_clode_root_dir",
    "initialize_runtime",
]