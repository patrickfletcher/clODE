from __future__ import annotations

from ..observers.types import ObserverParams
from ..problem.definition import ProblemInfo
from ..problem.source import RhsSource
from ..runtime.selection import RuntimeSelection
from ..simulation._protocols import FeatureBackend, SimulatorBackend, TrajectoryBackend
from .executors import (
    OpenCLFeatureExecutor,
    OpenCLTrajectoryExecutor,
    OpenCLTransientExecutor,
)
from .runtime import OpenCLRuntime


def _create_opencl_runtime(runtime_selection: RuntimeSelection | None) -> OpenCLRuntime:
    if runtime_selection is None:
        raise ValueError("OpenCL execution requires runtime selection metadata")

    device_id = runtime_selection.device_id
    if runtime_selection.device_ids is not None:
        if len(runtime_selection.device_ids) != 1:
            raise ValueError(
                "OpenCL execution currently supports exactly one selected device"
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
    clode_root: str,
    runtime_selection: RuntimeSelection | None = None,
) -> SimulatorBackend:
    return OpenCLTransientExecutor(
        problem_info,
        rhs_source,
        stepper,
        single_precision,
        _create_opencl_runtime(runtime_selection),
        clode_root,
    )


def create_trajectory_backend(
    problem_info: ProblemInfo,
    rhs_source: RhsSource,
    stepper: str,
    single_precision: bool,
    clode_root: str,
    runtime_selection: RuntimeSelection | None = None,
) -> TrajectoryBackend:
    return OpenCLTrajectoryExecutor(
        problem_info,
        rhs_source,
        stepper,
        single_precision,
        _create_opencl_runtime(runtime_selection),
        clode_root,
    )


def create_feature_backend(
    problem_info: ProblemInfo,
    rhs_source: RhsSource,
    stepper: str,
    observer: str,
    observer_params: ObserverParams,
    single_precision: bool,
    clode_root: str,
    runtime_selection: RuntimeSelection | None = None,
) -> FeatureBackend:
    return OpenCLFeatureExecutor(
        problem_info,
        rhs_source,
        stepper,
        observer,
        observer_params,
        single_precision,
        _create_opencl_runtime(runtime_selection),
        clode_root,
    )


__all__ = [
    "RuntimeSelection",
    "create_feature_backend",
    "create_simulator_backend",
    "create_trajectory_backend",
]