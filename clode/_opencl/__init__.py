from .buffers import ArrayLayout, BufferManager, CommonBuffers, FeatureBuffers, TrajectoryBuffers
from .errors import (
    BuildError,
    DoublePrecisionNotSupportedError,
    PyOpenCLBackendError,
    PyOpenCLDependencyError,
    PyOpenCLValidationError,
    RegistryValidationError,
    RhsValidationError,
    UnsupportedObserverError,
    UnsupportedStepperError,
)
from .executors import (
    PyOpenCLFeatureBackend,
    PyOpenCLTrajectoryBackend,
    PyOpenCLTransientBackend,
)
from .models import (
    PYOPENCL_BACKEND_VERSION,
    BuildKey,
    KernelKind,
    Precision,
    ProblemShape,
    ProgramBundle,
    RhsSource,
    SourceBundle,
)
from .registry import KernelRegistry
from .runtime import OpenCLRuntime
from .source_builder import SourceBuilder
from .observer_metadata import ObserverMetadata, get_observer_metadata, is_two_pass_observer
from .structs import (
    MatchedStruct,
    get_observer_params_struct,
    get_solver_params_struct,
    match_struct_dtype,
    pack_observer_params,
    pack_solver_params,
)
from .program_cache import ProgramCache

OPENCL_BACKEND_VERSION = PYOPENCL_BACKEND_VERSION
OpenCLFeatureExecutor = PyOpenCLFeatureBackend
OpenCLTrajectoryExecutor = PyOpenCLTrajectoryBackend
OpenCLTransientExecutor = PyOpenCLTransientBackend

__all__ = [
    "ArrayLayout",
    "BuildError",
    "BuildKey",
    "BufferManager",
    "CommonBuffers",
    "DoublePrecisionNotSupportedError",
    "FeatureBuffers",
    "MatchedStruct",
    "match_struct_dtype",
    "KernelKind",
    "KernelRegistry",
    "OPENCL_BACKEND_VERSION",
    "OpenCLRuntime",
    "OpenCLFeatureExecutor",
    "OpenCLTrajectoryExecutor",
    "OpenCLTransientExecutor",
    "ObserverMetadata",
    "PYOPENCL_BACKEND_VERSION",
    "Precision",
    "ProgramCache",
    "ProblemShape",
    "ProgramBundle",
    "PyOpenCLFeatureBackend",
    "PyOpenCLTrajectoryBackend",
    "PyOpenCLTransientBackend",
    "PyOpenCLBackendError",
    "PyOpenCLDependencyError",
    "PyOpenCLValidationError",
    "RegistryValidationError",
    "RhsSource",
    "RhsValidationError",
    "SourceBuilder",
    "SourceBundle",
    "TrajectoryBuffers",
    "get_observer_metadata",
    "get_observer_params_struct",
    "get_solver_params_struct",
    "is_two_pass_observer",
    "pack_observer_params",
    "pack_solver_params",
    "UnsupportedObserverError",
    "UnsupportedStepperError",
]