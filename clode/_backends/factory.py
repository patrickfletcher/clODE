from __future__ import annotations

from dataclasses import dataclass

from .._pyopencl.executors import (
    PyOpenCLFeatureBackend,
    PyOpenCLTrajectoryBackend,
    PyOpenCLTransientBackend,
)
from .._pyopencl.runtime import OpenCLRuntime
from ..problem.definition import ProblemInfo
from ..runtime import CLDeviceType, CLVendor, OpenCLResource, resolve_backend_name
from ..types import ObserverParams
from .protocol import FeatureBackend, SimulatorBackend, TrajectoryBackend
from .rhs import RhsSource


@dataclass(frozen=True, slots=True)
class RuntimeSelection:
    device_type: CLDeviceType | None
    vendor: CLVendor | None
    platform_id: int | None
    device_id: int | None
    device_ids: tuple[int, ...] | None
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
    resolve_backend_name(backend_name)
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
    resolve_backend_name(backend_name)
    return PyOpenCLTrajectoryBackend(
        problem_info,
        rhs_source,
        stepper,
        single_precision,
        _create_pyopencl_runtime(runtime_selection),
        clode_root,
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
    resolve_backend_name(backend_name)
    return PyOpenCLFeatureBackend(
        problem_info,
        rhs_source,
        stepper,
        observer,
        observer_params,
        single_precision,
        _create_pyopencl_runtime(runtime_selection),
        clode_root,
    )
