from __future__ import annotations

import os
from typing import Final

from clode.cpp.clode_cpp_wrapper import ObserverParams, ProblemInfo

from ..runtime import OpenCLResource
from .cpp import CppFeatureBackend, CppSimulatorBackend, CppTrajectoryBackend
from .protocol import FeatureBackend, SimulatorBackend, TrajectoryBackend

_DEFAULT_BACKEND: Final[str] = "cpp"
_BACKEND_ENVVAR: Final[str] = "_CLODE_BACKEND"


def _resolve_backend_name(backend_name: str | None = None) -> str:
    resolved = backend_name or os.getenv(_BACKEND_ENVVAR, _DEFAULT_BACKEND)
    if resolved != "cpp":
        raise ValueError(
            f"Unsupported clODE backend '{resolved}'. Supported backends: ['cpp']"
        )
    return resolved


def create_simulator_backend(
    problem_info: ProblemInfo,
    stepper: str,
    single_precision: bool,
    runtime: OpenCLResource,
    clode_root: str,
    backend_name: str | None = None,
) -> SimulatorBackend:
    _resolve_backend_name(backend_name)
    return CppSimulatorBackend(
        problem_info, stepper, single_precision, runtime, clode_root
    )


def create_trajectory_backend(
    problem_info: ProblemInfo,
    stepper: str,
    single_precision: bool,
    runtime: OpenCLResource,
    clode_root: str,
    backend_name: str | None = None,
) -> TrajectoryBackend:
    _resolve_backend_name(backend_name)
    return CppTrajectoryBackend(
        problem_info, stepper, single_precision, runtime, clode_root
    )


def create_feature_backend(
    problem_info: ProblemInfo,
    stepper: str,
    observer: str,
    observer_params: ObserverParams,
    single_precision: bool,
    runtime: OpenCLResource,
    clode_root: str,
    backend_name: str | None = None,
) -> FeatureBackend:
    _resolve_backend_name(backend_name)
    return CppFeatureBackend(
        problem_info,
        stepper,
        observer,
        observer_params,
        single_precision,
        runtime,
        clode_root,
    )
