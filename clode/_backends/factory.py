from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Final

from clode.cpp.clode_cpp_wrapper import ObserverParams, ProblemInfo

from .._pyopencl.executors import PyOpenCLTrajectoryBackend, PyOpenCLTransientBackend
from .._pyopencl.runtime import OpenCLRuntime
from ..runtime import CLDeviceType, CLVendor
from ..runtime import OpenCLResource
from .cpp import CppFeatureBackend, CppSimulatorBackend, CppTrajectoryBackend
from .protocol import FeatureBackend, SimulatorBackend, TrajectoryBackend
from .rhs import RhsSource

_DEFAULT_BACKEND: Final[str] = "cpp"
_BACKEND_ENVVAR: Final[str] = "_CLODE_BACKEND"


@dataclass(frozen=True, slots=True)
class RuntimeSelection:
    device_type: CLDeviceType | None
    vendor: CLVendor | None
    platform_id: int | None
    device_id: int | None
    device_ids: tuple[int, ...] | None


def _resolve_backend_name(backend_name: str | None = None) -> str:
    resolved = backend_name or os.getenv(_BACKEND_ENVVAR, _DEFAULT_BACKEND)
    if resolved not in {"cpp", "pyopencl"}:
        raise ValueError(
            f"Unsupported clODE backend '{resolved}'. Supported backends: ['cpp', 'pyopencl']"
        )
    return resolved


def _create_pyopencl_runtime(runtime_selection: RuntimeSelection | None) -> OpenCLRuntime:
    if runtime_selection is None:
        raise ValueError("PyOpenCL backend requires runtime selection metadata")

    device_id = runtime_selection.device_id
    if runtime_selection.device_ids is not None:
        if len(runtime_selection.device_ids) != 1:
            raise ValueError(
                "PyOpenCL backend currently supports exactly one selected device"
            )
        device_id = runtime_selection.device_ids[0]

    return OpenCLRuntime.create(
        device_type=runtime_selection.device_type,
        vendor=runtime_selection.vendor,
        platform_id=runtime_selection.platform_id,
        device_id=device_id,
    )


def create_simulator_backend(
    problem_info: ProblemInfo,
    rhs_source: RhsSource,
    stepper: str,
    single_precision: bool,
    runtime: OpenCLResource,
    clode_root: str,
    runtime_selection: RuntimeSelection | None = None,
    backend_name: str | None = None,
) -> SimulatorBackend:
    resolved = _resolve_backend_name(backend_name)
    if resolved == "cpp":
        return CppSimulatorBackend(
            problem_info, rhs_source, stepper, single_precision, runtime, clode_root
        )
    return PyOpenCLTransientBackend(
        problem_info,
        rhs_source,
        stepper,
        single_precision,
        _create_pyopencl_runtime(runtime_selection),
        clode_root,
    )


def create_trajectory_backend(
    problem_info: ProblemInfo,
    rhs_source: RhsSource,
    stepper: str,
    single_precision: bool,
    runtime: OpenCLResource,
    clode_root: str,
    runtime_selection: RuntimeSelection | None = None,
    backend_name: str | None = None,
) -> TrajectoryBackend:
    resolved = _resolve_backend_name(backend_name)
    if resolved == "pyopencl":
        return PyOpenCLTrajectoryBackend(
            problem_info,
            rhs_source,
            stepper,
            single_precision,
            _create_pyopencl_runtime(runtime_selection),
            clode_root,
        )
    return CppTrajectoryBackend(
        problem_info, rhs_source, stepper, single_precision, runtime, clode_root
    )


def create_feature_backend(
    problem_info: ProblemInfo,
    rhs_source: RhsSource,
    stepper: str,
    observer: str,
    observer_params: ObserverParams,
    single_precision: bool,
    runtime: OpenCLResource,
    clode_root: str,
    runtime_selection: RuntimeSelection | None = None,
    backend_name: str | None = None,
) -> FeatureBackend:
    resolved = _resolve_backend_name(backend_name)
    if resolved == "pyopencl":
        raise NotImplementedError("PyOpenCL feature backend is not implemented yet")
    return CppFeatureBackend(
        problem_info,
        rhs_source,
        stepper,
        observer,
        observer_params,
        single_precision,
        runtime,
        clode_root,
    )
