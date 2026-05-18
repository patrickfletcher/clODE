from __future__ import annotations

import argparse
import os
import sys
import traceback
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Callable

REPO_ROOT = Path(__file__).resolve().parent.parent
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

import numpy as np
import pyopencl as cl
import pyopencl._cl as cl_module
import pyopencl.tools as cl_tools

from clode._opencl.models import Precision, ProblemShape
from clode._opencl.observer_metadata import (
    _observer_data_base_dtype,
    _observer_data_struct_name,
)
from clode.problem._core import ProblemInfo

ENV_VARS = (
    "CLODE_TEST_PLATFORM_ID",
    "CLODE_TEST_DEVICE_ID",
    "POCL_KERNELLIB_NAME",
    "POCL_LLVM_CPU_NAME",
    "POCL_KERNEL_CACHE",
    "PYOPENCL_COMPILER_OUTPUT",
)


def _env_int(name: str, fallback: int) -> int:
    raw_value = os.environ.get(name)
    if raw_value is None:
        return fallback
    try:
        return int(raw_value)
    except ValueError as exc:
        raise SystemExit(f"Environment variable {name} must be an integer, got {raw_value!r}") from exc


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Probe the selected OpenCL runtime with increasingly specific build/layout checks. "
            "This is useful for separating runtime/compiler failures from clODE kernel issues."
        )
    )
    parser.add_argument(
        "--platform-id",
        type=int,
        default=_env_int("CLODE_TEST_PLATFORM_ID", 0),
        help="Platform index to probe. Defaults to CLODE_TEST_PLATFORM_ID or 0.",
    )
    parser.add_argument(
        "--device-id",
        type=int,
        default=_env_int("CLODE_TEST_DEVICE_ID", 0),
        help="Device index within the selected platform. Defaults to CLODE_TEST_DEVICE_ID or 0.",
    )
    parser.add_argument(
        "--max-event-timestamps",
        type=int,
        default=3,
        help="Event-storage count used for the clODE localmax observer struct probe.",
    )
    parser.add_argument(
        "--list",
        action="store_true",
        help="List visible platforms and devices, then exit.",
    )
    parser.add_argument(
        "--traceback",
        action="store_true",
        help="Print Python tracebacks for failing probes.",
    )
    return parser


def _distribution_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "<not installed>"


def _print_runtime_metadata() -> None:
    print("Runtime metadata:")
    print(f"  pyopencl={_distribution_version('pyopencl')}")
    print(f"  pyopencl._cl={cl_module.__file__}")


def _device_type_name(device: cl.Device) -> str:
    try:
        return cl.device_type.to_string(device.type)
    except Exception:
        return str(device.type)


def _print_runtime_inventory() -> list[cl.Platform]:
    platforms = cl.get_platforms()
    if not platforms:
        raise SystemExit("No OpenCL platforms were found.")

    print("Visible OpenCL platforms and devices:")
    for platform_index, platform in enumerate(platforms):
        print(
            f"- Platform {platform_index}: {platform.name} | vendor={platform.vendor} | version={platform.version}"
        )
        devices = platform.get_devices()
        if not devices:
            print("  (no devices)")
            continue
        for device_index, device in enumerate(devices):
            print(
                "  "
                f"Device {device_index}: {device.name} | type={_device_type_name(device)} | "
                f"version={device.version} | driver={device.driver_version}"
            )
    return platforms


def _select_device(platform_id: int, device_id: int) -> tuple[cl.Platform, cl.Device]:
    platforms = cl.get_platforms()
    if platform_id < 0 or platform_id >= len(platforms):
        raise SystemExit(
            f"platform-id {platform_id} is out of range; found {len(platforms)} platform(s)."
        )

    platform = platforms[platform_id]
    devices = platform.get_devices()
    if not devices:
        raise SystemExit(f"Platform {platform_id} exposes no devices.")
    if device_id < 0 or device_id >= len(devices):
        raise SystemExit(
            f"device-id {device_id} is out of range for platform {platform_id}; "
            f"found {len(devices)} device(s)."
        )
    return platform, devices[device_id]


def _print_selected_device(platform_id: int, device_id: int, platform: cl.Platform, device: cl.Device) -> None:
    print("")
    print(f"Selected platform/device: {platform_id}/{device_id}")
    print(f"  platform: {platform.name}")
    print(f"  vendor: {platform.vendor}")
    print(f"  platform version: {platform.version}")
    print(f"  device: {device.name}")
    print(f"  device type: {_device_type_name(device)}")
    print(f"  device version: {device.version}")
    print(f"  driver version: {device.driver_version}")
    opencl_c_version = getattr(device, "opencl_c_version", None)
    if opencl_c_version:
        print(f"  OpenCL C version: {opencl_c_version}")


def _print_relevant_env() -> None:
    print("")
    print("Relevant environment variables:")
    for name in ENV_VARS:
        value = os.environ.get(name)
        rendered = "<unset>" if value is None else value
        print(f"  {name}={rendered}")


def _probe_trivial_kernel(context: cl.Context, device: cl.Device) -> None:
    cl.Program(context, "__kernel void noop(void) {}\n").build(devices=[device])


def _probe_diagnostic_struct(context: cl.Context, device: cl.Device) -> None:
    diagnostic_dtype = np.dtype(
        [("count", np.uint32), ("value", np.float64), ("flag", np.uint32)],
        align=True,
    )
    matched_dtype, c_decl = cl_tools.match_dtype_to_c_struct(
        device,
        "diagnostic_struct",
        diagnostic_dtype,
        context,
    )
    print(f"  struct name: diagnostic_struct")
    print(f"  numpy itemsize: {diagnostic_dtype.itemsize}")
    print(f"  matched itemsize: {matched_dtype.itemsize}")
    print("  c declaration:")
    for line in c_decl.strip().splitlines():
        print(f"    {line}")


def _probe_clode_localmax_struct(
    context: cl.Context,
    device: cl.Device,
    max_event_timestamps: int,
) -> None:
    problem_info = ProblemInfo(
        src_file="diagnostic.cl",
        vars=["x"],
        pars=["p"],
        aux=[],
        num_noise=0,
    )
    shape = ProblemShape.from_problem_info(problem_info)
    base_dtype = _observer_data_base_dtype(
        "localmax",
        shape,
        Precision.DOUBLE,
        max_event_timestamps,
    )
    struct_name = _observer_data_struct_name(
        "localmax",
        shape,
        Precision.DOUBLE,
        max_event_timestamps,
    )
    matched_dtype, _c_decl = cl_tools.match_dtype_to_c_struct(
        device,
        struct_name,
        base_dtype,
        context,
    )
    field_count = len(base_dtype.names or ())
    print(f"  struct name: {struct_name}")
    print(f"  fields: {field_count}")
    print(f"  numpy itemsize: {base_dtype.itemsize}")
    print(f"  matched itemsize: {matched_dtype.itemsize}")
    print(f"  max_event_timestamps: {max_event_timestamps}")


def _run_probe(
    title: str,
    probe: Callable[[], None],
    show_traceback: bool,
) -> bool:
    print("")
    print(f"== {title} ==")
    try:
        probe()
    except Exception as exc:
        print("status: FAIL")
        print(f"exception type: {exc.__class__.__name__}")
        print(f"exception message: {exc}")
        if show_traceback:
            traceback.print_exc(file=sys.stdout)
        return False

    print("status: OK")
    return True


def main(argv: list[str]) -> int:
    args = _build_parser().parse_args(argv)
    _print_runtime_metadata()
    print("")
    platforms = _print_runtime_inventory()
    if args.list:
        return 0
    if not platforms:
        return 1

    platform, device = _select_device(args.platform_id, args.device_id)
    _print_selected_device(args.platform_id, args.device_id, platform, device)
    _print_relevant_env()

    try:
        context = cl.Context([device])
    except Exception as exc:
        print("")
        print("Failed to create a context for the selected device.")
        print(f"exception type: {exc.__class__.__name__}")
        print(f"exception message: {exc}")
        if args.traceback:
            traceback.print_exc(file=sys.stdout)
        return 1

    probe_results = [
        _run_probe(
            "Probe 1: trivial kernel build",
            lambda: _probe_trivial_kernel(context, device),
            args.traceback,
        ),
        _run_probe(
            "Probe 2: PyOpenCL diagnostic struct layout",
            lambda: _probe_diagnostic_struct(context, device),
            args.traceback,
        ),
        _run_probe(
            "Probe 3: clODE localmax observer struct layout",
            lambda: _probe_clode_localmax_struct(
                context,
                device,
                args.max_event_timestamps,
            ),
            args.traceback,
        ),
    ]

    print("")
    if all(probe_results):
        print("All probes succeeded.")
        return 0

    failed_count = sum(1 for success in probe_results if not success)
    print(f"{failed_count} probe(s) failed.")
    return 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
