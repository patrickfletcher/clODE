# PyOpenCL Backend Design And Implementation Plan

## Status

- Status: In progress
- Audience: maintainers and contributors implementing the PyOpenCL migration
- Related document: `tmp/cpp_opencl_layer_audit.md`
- Completed groundwork: scope lock, authoritative migration test gate, internal backend seam, C++ backend adapter, and RHS source object integration

## Decision Summary

clODE should migrate to a Python-owned OpenCL runtime built on PyOpenCL while preserving the current public Python API and preserving the existing OpenCL numerical kernels as the behavioral reference.

The migration should not begin by unifying `transient`, `trajectory`, and `features` into a single new abstraction inside the current C++ layer. That is a valid long-term direction, but it is not the right first move. The first move is to:

1. introduce an internal backend boundary in Python
2. keep the current public API stable
3. implement a PyOpenCL backend behind that boundary
4. land transient parity first
5. add trajectory parity second
6. add feature parity third
7. only then simplify duplicated execution paths and revisit the broader monitor or sink abstraction

This sequencing minimizes risk, avoids redesigning code that is about to be replaced, and creates a clear path to shipping incremental value.

## Problem Statement

The current C++ layer successfully runs the numerical solver, but it is carrying concerns that are not unique to clODE:

- OpenCL context and queue management
- program source assembly
- build option generation
- program rebuild invalidation
- buffer allocation and data transfer
- pybind11 adapter glue

Those responsibilities duplicate logic already present in the Python wrappers and are the main source of maintenance cost. The kernels themselves are the part worth preserving.

At the same time, the current system relies on compile-time specialization for:

- precision
- stepper
- observer
- problem dimensions
- stored event count

That specialization is deeply embedded in the kernels and should be preserved during the first migration phase rather than removed immediately.

## Goals

### Primary goals

- Replace the C++ runtime layer with a Python runtime built on PyOpenCL.
- Preserve the current public Python API surface for `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` during the migration.
- Preserve the current OpenCL kernels and current behavioral semantics during the parity phases.
- Centralize build-key computation, source assembly, buffer allocation, and rebuild policy in Python.
- Make the system testable with side-by-side C++ backend vs PyOpenCL backend parity checks.

### Secondary goals

- Make the generated OpenCL program source easy to inspect and debug.
- Remove mixed host and kernel metadata concerns from the OpenCL source tree over time.
- Set up the codebase so later simplifications, including a broader monitor or sink abstraction, are straightforward.

## Non-goals

These items are explicitly out of scope for the first shipping version of the PyOpenCL backend.

- Redesigning the public Python constructor signatures or return types.
- Replacing compile-time dimension specialization with runtime-generic kernels.
- Real multi-device execution.
- A full unification of transient, trajectory, and feature execution into one kernel contract.
- New numerical methods, new observers, or changes to output semantics.
- Fixing every historical API inconsistency in the same change set.

## Preservation Rules

The migration should preserve the following contracts until parity is established and proven.

- `Simulator`, `TrajectorySimulator`, and `FeatureSimulator` remain the public entry points.
- `Stepper` and `Observer` enums remain stable.
- `TrajectoryOutput` and `ObserverOutput` remain stable in shape and semantics.
- Constructor defaults remain stable, even if they are not ideal.
- `get_program_string()` remains available.
- `seed_rng`, `shift_x0`, `shift_tspan`, `get_dt`, and `get_final_time` remain available.
- Single vs double precision behavior remains compile-time specialized.
- The current Fortran-order flattening and reshape semantics remain the backend implementation default.

## Design Principles

1. Preserve behavior first, simplify second.
2. Move ownership of runtime concerns to Python early.
3. Minimize the number of new abstractions introduced per phase.
4. Make rebuild logic explicit and deterministic.
5. Keep the public API stable while the backend evolves.
6. Prefer side-by-side parity tests over assumptions.
7. Defer conceptual cleanup until the PyOpenCL backend is already shipping.

## Proposed Target Architecture

### High-level shape

The public API stays where it is. The implementation underneath it changes from a direct dependency on the pybind C++ classes to an internal backend abstraction.

The architecture is:

```text
Public API
  Simulator / TrajectorySimulator / FeatureSimulator
      |
      v
Internal Backend Boundary
  SimulatorBackend / TrajectoryBackend / FeatureBackend
      |
      +--> Cpp backend adapter (temporary, for parity and rollout)
      |
      +--> PyOpenCL backend
              |
              +--> Runtime
              +--> Kernel registry
              +--> Source builder
              +--> Program cache
              +--> Buffer manager
              +--> Executor classes
```

Current state:

- the public simulators now delegate backend construction through `clode/_backends/factory.py`
- the active implementation path is the C++ adapter in `clode/_backends/cpp.py`
- public simulators now prepare an internal `RhsSource` object before backend construction
- the next backend milestone is PR 5: PyOpenCL core models and errors

### Recommended module layout

The public modules remain in place. New implementation modules are introduced under internal packages.

```text
clode/
  runtime.py
  solver.py
  trajectory.py
  features.py
  function_converter.py
  xpp_parser.py

  _backends/
    __init__.py
    protocol.py
    factory.py
    cpp.py

  _pyopencl/
    __init__.py
    models.py
    runtime.py
    registry.py
    source_builder.py
    program_cache.py
    buffers.py
    executors.py
    errors.py
```

This layout is intentionally conservative. It introduces enough separation to keep the implementation understandable without exploding the number of modules.

## Module Boundaries

### Public modules

#### `clode/runtime.py`

Responsibility:

- keep the public runtime API stable
- expose `initialize_runtime`, `print_opencl`, device info types, and logging helpers
- stop being the place where solver execution is implemented

Phase behavior:

- during transition, it may still expose the current C++ `OpenCLResource`
- once the PyOpenCL runtime is ready, it should expose a Python `OpenCLResource` compatibility class with the same public behavior

#### `clode/solver.py`, `clode/trajectory.py`, `clode/features.py`

Responsibility:

- remain the user-facing orchestration layer
- hold argument validation, ensemble shaping, and output reshaping
- delegate all device and kernel execution work to the backend object

They should no longer instantiate pybind classes directly once the backend boundary is in place.

### Internal backend boundary

#### `clode/_backends/protocol.py`

Responsibility:

- define the internal protocol used by the public simulators
- make backend switching possible without changing the public simulators again later

#### `clode/_backends/factory.py`

Responsibility:

- choose which backend implementation to instantiate
- allow an internal-only backend override for parity testing and staged rollout

#### `clode/_backends/cpp.py`

Responsibility:

- wrap the existing pybind-backed runtime behind the new protocol
- provide a stable reference backend during migration

This adapter is not the long-term design target. It is the migration stabilizer.

### PyOpenCL implementation modules

#### `clode/_pyopencl/models.py`

Responsibility:

- define small immutable data models and enums used throughout the PyOpenCL backend
- avoid spreading ad hoc tuples and string keys across the implementation

Expected contents:

- `Precision`
- `KernelKind`
- `ProblemShape`
- `RhsSource`
- `BuildKey`
- `SourceBundle`
- `ProgramBundle`

#### `clode/_pyopencl/runtime.py`

Responsibility:

- context creation
- queue creation
- device selection
- device capability reporting
- runtime-scoped caches

This is the Python replacement for the OpenCL parts of `OpenCLResource`.

#### `clode/_pyopencl/registry.py`

Responsibility:

- stepper registry
- observer registry
- observer metadata builders
- entrypoint-to-source mapping

It moves host metadata out of the `.cl` files over time. In phase one it may still duplicate some metadata from the current source tree because that is cheaper and safer than trying to introspect it from OpenCL files.

#### `clode/_pyopencl/source_builder.py`

Responsibility:

- assemble the program source for a requested build
- validate the RHS source
- compute a stable source digest
- return both source text and build options in one object

This module should own the canonical source assembly path.

#### `clode/_pyopencl/program_cache.py`

Responsibility:

- compile PyOpenCL programs on demand
- cache compiled programs by build key
- cache kernel handles associated with a built program
- capture and surface build errors with the source text and build options attached

#### `clode/_pyopencl/buffers.py`

Responsibility:

- allocate and reuse device buffers
- flatten and reshape arrays consistently
- own the implementation of the current Fortran-order layout contract
- keep host-to-device and device-to-host transfer code out of the executor classes

#### `clode/_pyopencl/executors.py`

Responsibility:

- implement the three backend executors that match the current public simulator behaviors
- share setup logic where it is actually common
- keep execution-specific buffer sets isolated

The important design choice is that these executors stay separate through the parity phases.

## Backend Boundary Design

The backend boundary should mirror the operations already assumed by the public Python classes.

### Protocols

```python
from __future__ import annotations

from typing import Protocol, Optional


class SimulatorBackend(Protocol):
    def build(self) -> None: ...
    def set_problem_data(self, x0: list[float], pars: list[float]) -> None: ...
    def set_x0(self, x0: list[float]) -> None: ...
    def set_pars(self, pars: list[float]) -> None: ...
    def set_tspan(self, tspan: tuple[float, float]) -> None: ...
    def set_solver_params(self, params: object) -> None: ...
    def seed_rng(self, seed: Optional[int] = None) -> None: ...
    def transient(self) -> None: ...
    def shift_x0(self) -> None: ...
    def shift_tspan(self) -> None: ...
    def get_x0(self) -> list[float]: ...
    def get_xf(self) -> list[float]: ...
    def get_dt(self) -> list[float]: ...
    def get_tf(self) -> list[float]: ...
    def get_tspan(self) -> tuple[float, float]: ...
    def get_program_string(self) -> str: ...
    def get_available_steppers(self) -> list[str]: ...
    def print_status(self) -> None: ...


class TrajectoryBackend(SimulatorBackend, Protocol):
    def trajectory(self) -> None: ...
    def get_t(self) -> list[float]: ...
    def get_x(self) -> list[float]: ...
    def get_dx(self) -> list[float]: ...
    def get_aux(self) -> list[float]: ...
    def get_n_stored(self) -> list[int]: ...


class FeatureBackend(SimulatorBackend, Protocol):
    def set_observer(self, observer: str) -> None: ...
    def set_observer_params(self, params: object) -> None: ...
    def initialize_observer(self) -> None: ...
    def is_observer_initialized(self) -> bool: ...
    def features(self, reinitialize_observer: Optional[bool] = None) -> None: ...
    def get_f(self) -> list[float]: ...
    def get_n_features(self) -> int: ...
    def get_feature_names(self) -> list[str]: ...
    def get_available_observers(self) -> list[str]: ...
```

The public simulator classes already assume something close to this shape. Making it explicit is the key enabling move.

### C++ adapter skeleton

```python
class CppSimulatorBackend:
    def __init__(self, problem_info, stepper, single_precision, runtime, kernel_root):
        self._integrator = SimulatorBase(
            problem_info,
            stepper,
            single_precision,
            runtime,
            kernel_root,
        )

    def build(self) -> None:
        self._integrator.build_cl()

    def set_problem_data(self, x0: list[float], pars: list[float]) -> None:
        self._integrator.set_problem_data(x0, pars)

    # remaining methods delegate directly
```

This adapter should be small and boring.

## PyOpenCL Core Models

These are the minimum data models needed to avoid key logic being encoded as positional tuples and free-form dictionaries.

```python
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import Optional


class Precision(Enum):
    SINGLE = "single"
    DOUBLE = "double"


class KernelKind(Enum):
    TRANSIENT = "transient"
    TRAJECTORY = "trajectory"
    FEATURES = "features"


@dataclass(frozen=True)
class ProblemShape:
    n_var: int
    n_par: int
    n_aux: int
    n_wiener: int


@dataclass(frozen=True)
class RhsSource:
    origin_label: str
    text: str
    digest: str


@dataclass(frozen=True)
class BuildKey:
    backend_version: str
    kernel_kind: KernelKind
    precision: Precision
    stepper_name: str
    observer_name: Optional[str]
    problem_shape: ProblemShape
    n_store_events: int
    rhs_digest: str
    kernel_tree_digest: str
    debug_build: bool = False


@dataclass(frozen=True)
class SourceBundle:
    build_key: BuildKey
    source_text: str
    build_options: tuple[str, ...]
    kernel_names: tuple[str, ...]


@dataclass
class ProgramBundle:
    build_key: BuildKey
    source_bundle: SourceBundle
    program: object
    kernels: dict[str, object]
```

### Why the build key includes digests

The build key must include more than dimensions and stepper selection.

It also needs:

- an RHS digest, because the program changes when the problem source changes
- a kernel tree digest, because internal kernel library edits must invalidate the cache even if the user problem did not change
- a backend version, because the Python source builder itself may change output semantics over time

This is necessary to make cache hits correct.

## Runtime Design

The runtime object owns the OpenCL context and queue and carries a program cache scoped to that context.

### Runtime skeleton

```python
class OpenCLRuntime:
    def __init__(self, context, queue, platform_info, device_info):
        self.context = context
        self.queue = queue
        self.platform_info = platform_info
        self.device_info = device_info
        self.program_cache = ProgramCache()

    @classmethod
    def create(
        cls,
        device_type=None,
        vendor=None,
        platform_id=None,
        device_id=None,
    ) -> "OpenCLRuntime":
        # explicit single-device selection in phase one
        ...

    def get_double_support(self) -> bool:
        ...

    def get_max_memory_alloc_size(self) -> int:
        ...

    def get_device_cl_version(self) -> str:
        ...
```

### Compatibility wrapper

The public `initialize_runtime(...)` function should continue to return an `OpenCLResource`-like object. During rollout there are two valid approaches.

Approach A:

- keep returning the current C++ object while public simulators still default to the C++ backend
- use the Python runtime internally only when explicitly selecting the PyOpenCL backend

Approach B:

- introduce a Python `OpenCLResource` compatibility class early
- make both backends able to work with it through adapters

Approach A is lower risk for the initial migration.

## Source Builder Design

The source builder must be deterministic, inspectable, and phaseable.

### Source builder responsibilities

- read the entrypoint kernel files for the requested kernel kind
- read the RHS source
- compute the kernel tree digest
- compute the build key
- assemble the final source text
- emit the build options tuple
- validate that the RHS defines `getRHS(...)`
- provide the exact program text exposed by `get_program_string()`

### Phase-one source strategy

Phase one should preserve the current source organization.

- `transient` builds from `transient.cl + rhs`
- `trajectory` builds from `transient.cl + trajectory.cl + rhs`
- `features` builds from `transient.cl + initializeObserver.cl + features.cl + rhs`

The builder should continue to rely on `#include` directives for the shared kernel library and should pass the kernel root as an include directory.

This is the minimal-change path.

### Phase-two source strategy

After parity, the builder may optionally resolve local includes itself and emit a fully expanded monolithic source string. That can remove dependence on compiler include behavior and make metadata extraction easier, but it should not be the first migration step.

### Source builder skeleton

```python
class SourceBuilder:
    def __init__(self, kernel_root: Path, registry: "KernelRegistry"):
        self._kernel_root = kernel_root
        self._registry = registry

    def build(
        self,
        kernel_kind: KernelKind,
        precision: Precision,
        stepper_name: str,
        problem_shape: ProblemShape,
        rhs: RhsSource,
        observer_name: str | None = None,
        n_store_events: int = 0,
        debug_build: bool = False,
    ) -> SourceBundle:
        self._registry.validate_stepper(stepper_name)
        self._registry.validate_observer(observer_name)
        self._validate_rhs(rhs)

        entrypoint_paths = self._registry.get_entrypoint_paths(kernel_kind)
        entrypoint_source = "".join(path.read_text() for path in entrypoint_paths)
        kernel_tree_digest = self._compute_kernel_tree_digest(entrypoint_paths)

        build_key = BuildKey(
            backend_version="1",
            kernel_kind=kernel_kind,
            precision=precision,
            stepper_name=stepper_name,
            observer_name=observer_name,
            problem_shape=problem_shape,
            n_store_events=n_store_events,
            rhs_digest=rhs.digest,
            kernel_tree_digest=kernel_tree_digest,
            debug_build=debug_build,
        )

        build_options = self._make_build_options(build_key)
        source_text = entrypoint_source + rhs.text
        kernel_names = self._registry.get_kernel_names(kernel_kind)

        return SourceBundle(
            build_key=build_key,
            source_text=source_text,
            build_options=build_options,
            kernel_names=kernel_names,
        )

    def _make_build_options(self, key: BuildKey) -> tuple[str, ...]:
        options = []
        options.append("-DCLODE_SINGLE_PRECISION" if key.precision is Precision.SINGLE else "-DCLODE_DOUBLE_PRECISION")
        options.append(f"-D{self._registry.get_stepper_define(key.stepper_name)}")
        options.append(f"-DN_VAR={key.problem_shape.n_var}")
        options.append(f"-DN_PAR={key.problem_shape.n_par}")
        options.append(f"-DN_AUX={key.problem_shape.n_aux}")
        options.append(f"-DN_WIENER={key.problem_shape.n_wiener}")
        options.append(f"-I{self._kernel_root}")
        if key.observer_name is not None:
            options.append(f"-D{self._registry.get_observer_define(key.observer_name)}")
            options.append(f"-DN_STORE_EVENTS={key.n_store_events}")
        return tuple(options)

    def _validate_rhs(self, rhs: RhsSource) -> None:
        # phase-one validation can be regex-based and conservative
        ...
```

### Registry skeleton

```python
class KernelRegistry:
    _STEPPERS = {
        "euler": "EXPLICIT_EULER",
        "heun": "EXPLICIT_HEUN",
        "rk4": "EXPLICIT_RK4",
        "bs23": "EXPLICIT_BS23",
        "dopri5": "EXPLICIT_DOPRI5",
        "seuler": "STOCHASTIC_EULER",
    }

    _OBSERVERS = {
        "basic": "USE_OBSERVER_BASIC",
        "basicall": "USE_OBSERVER_BASIC_ALLVAR",
        "localmax": "USE_OBSERVER_LOCAL_MAX",
        "nhood1": "USE_OBSERVER_NHOOD_1",
        "nhood2": "USE_OBSERVER_NHOOD_2",
        "thresh2": "USE_OBSERVER_THRESHOLD_2",
    }

    _ENTRYPOINTS = {
        KernelKind.TRANSIENT: ("transient.cl",),
        KernelKind.TRAJECTORY: ("transient.cl", "trajectory.cl"),
        KernelKind.FEATURES: ("transient.cl", "initializeObserver.cl", "features.cl"),
    }

    _KERNEL_NAMES = {
        KernelKind.TRANSIENT: ("transient",),
        KernelKind.TRAJECTORY: ("transient", "trajectory"),
        KernelKind.FEATURES: ("transient", "initializeObserver", "features"),
    }
```

The registry may begin as static Python data. That is acceptable in phase one.

## Program Cache Design

The program cache is runtime-scoped. A compiled `pyopencl.Program` is valid only for the context and devices it was built against. Therefore, the simplest correct design is:

- one cache per `OpenCLRuntime`
- cache keyed by `BuildKey`
- no attempt at cross-runtime compiled-binary reuse in phase one

### Program cache responsibilities

- perform cache lookup by build key
- compile on cache miss
- create and cache requested kernels
- surface build errors with source and options attached
- optionally record build timings for diagnostics

### Program cache skeleton

```python
class BuildError(RuntimeError):
    def __init__(self, message: str, source_text: str, build_options: tuple[str, ...], build_log: str):
        super().__init__(message)
        self.source_text = source_text
        self.build_options = build_options
        self.build_log = build_log


class ProgramCache:
    def __init__(self):
        self._cache: dict[BuildKey, ProgramBundle] = {}

    def get_or_build(self, runtime: OpenCLRuntime, source_bundle: SourceBundle) -> ProgramBundle:
        existing = self._cache.get(source_bundle.build_key)
        if existing is not None:
            return existing

        try:
            program = cl.Program(runtime.context, source_bundle.source_text).build(
                options=list(source_bundle.build_options)
            )
        except Exception as exc:
            build_log = self._extract_build_log(exc)
            raise BuildError(
                "OpenCL program build failed",
                source_text=source_bundle.source_text,
                build_options=source_bundle.build_options,
                build_log=build_log,
            ) from exc

        kernels = {
            kernel_name: cl.Kernel(program, kernel_name)
            for kernel_name in source_bundle.kernel_names
        }
        bundle = ProgramBundle(
            build_key=source_bundle.build_key,
            source_bundle=source_bundle,
            program=program,
            kernels=kernels,
        )
        self._cache[source_bundle.build_key] = bundle
        return bundle

    def clear(self) -> None:
        self._cache.clear()

    def _extract_build_log(self, exc: Exception) -> str:
        ...
```

### Cache invalidation policy

Cache invalidation happens structurally, not procedurally.

If any field that affects generated code changes, the build key changes, so the cache misses naturally.

This removes the need for partially maintained `_cl_program_is_valid` flags inside public wrappers.

## Buffer Manager Design

The buffer manager is responsible for the device-side state model and the reshape rules currently split across C++ and Python.

### Responsibilities

- own the current flatten order
- allocate common buffers based on ensemble size and precision
- allocate optional trajectory buffers
- allocate optional feature buffers
- provide typed upload and download helpers
- preserve the existing public shape semantics for returned NumPy arrays

### Design choice

The buffer manager should not try to be generic over every possible kernel shape. It should explicitly support the current solver families.

Recommended objects:

- `CommonBuffers`
- `TrajectoryBuffers`
- `FeatureBuffers`
- `ArrayLayout`

### Skeleton

```python
@dataclass
class CommonBuffers:
    tspan: object
    solver_params: object
    x0: object
    pars: object
    xf: object
    rng_state: object
    dt: object
    tf: object


class ArrayLayout:
    @staticmethod
    def flatten_problem_matrix(array: np.ndarray) -> np.ndarray:
        return np.asarray(array, dtype=np.float64).flatten(order="F")

    @staticmethod
    def reshape_state(array: list[float], ensemble_size: int, width: int) -> np.ndarray:
        return np.asarray(array, dtype=np.float64).reshape((ensemble_size, width), order="F")


class BufferManager:
    def __init__(self, runtime: OpenCLRuntime, precision: Precision):
        self._runtime = runtime
        self._precision = precision

    def allocate_common(self, ensemble_size: int, shape: ProblemShape) -> CommonBuffers:
        ...

    def allocate_trajectory(self, ensemble_size: int, shape: ProblemShape, max_store: int) -> object:
        ...

    def allocate_features(self, ensemble_size: int, n_features: int, observer_data_size: int) -> object:
        ...
```

## Executor Design

The executor layer is where behavior parity lives. This is not the place to impose the future monitor or sink abstraction yet.

### Why keep three executors in the parity phases

Although `transient`, `trajectory`, and `features` are conceptually related, they differ in all of the following today.

- kernel entrypoints
- build keys
- output buffers
- optional observer warmup
- persistent observer state
- output reshaping

Trying to force them into one executor before parity would recreate the failed or unfinished unification work already visible in `odedriver.cl`.

### Recommended executors

- `PyOpenCLTransientBackend`
- `PyOpenCLTrajectoryBackend`
- `PyOpenCLFeatureBackend`

These classes should share small reusable helpers for common setup but remain separate classes.

### Base executor skeleton

```python
class PyOpenCLSimulatorBackend:
    def __init__(self, problem_info, stepper, single_precision, runtime, kernel_root, rhs_source):
        self._problem_info = problem_info
        self._stepper = stepper
        self._precision = Precision.SINGLE if single_precision else Precision.DOUBLE
        self._runtime = runtime
        self._kernel_root = kernel_root
        self._rhs_source = rhs_source
        self._builder = SourceBuilder(kernel_root, KernelRegistry())
        self._buffers = BufferManager(runtime, self._precision)
        self._program_bundle = None

    def build(self) -> None:
        source_bundle = self._builder.build(
            kernel_kind=KernelKind.TRANSIENT,
            precision=self._precision,
            stepper_name=self._stepper,
            problem_shape=self._problem_shape(),
            rhs=self._rhs_source,
        )
        self._program_bundle = self._runtime.program_cache.get_or_build(self._runtime, source_bundle)

    def _problem_shape(self) -> ProblemShape:
        return ProblemShape(
            n_var=self._problem_info.num_var,
            n_par=self._problem_info.num_par,
            n_aux=self._problem_info.num_aux,
            n_wiener=self._problem_info.num_noise,
        )
```

Trajectory and feature backends then extend this with their additional buffers and kernel launches.

## Public API Integration Plan

The public simulators should gradually stop depending on pybind implementation details.

### Required refactor

`solver.py`, `trajectory.py`, and `features.py` should call a backend factory instead of directly instantiating pybind types.

For example:

```python
from ._backends.factory import create_simulator_backend


def _create_integrator(self) -> None:
    self._integrator = create_simulator_backend(
        kind="simulator",
        problem_info=self._pi,
        stepper=self._stepper.value,
        single_precision=self._single_precision,
        runtime=self._runtime,
        kernel_root=_clode_root_dir,
        rhs_source=self._rhs_source,
    )
```

This is the most important structural refactor in the plan. It is also small enough to land before any PyOpenCL kernel execution is implemented.

### RHS source preparation

The public solver layer should stop treating the RHS as only a file path. Internally it should produce an `RhsSource` object that always contains source text and a digest.

Public behavior can remain file-based in phase one. Internally the backend should always see:

- source text
- an origin label
- a stable digest

This change makes caching and build diagnostics much cleaner.

## Error Handling And Diagnostics

The PyOpenCL backend should be stricter and more informative than the current backend.

### Minimum diagnostics to add

- explicit exception type for program build failures
- build log included on failure
- source text included on failure
- build options included on failure
- validation error when RHS does not define `getRHS(...)`
- validation error for unsupported stepper or observer names

### `get_program_string()` behavior

The backend should return the exact source string handed to the compiler plus a formatted view of build options. That matches the spirit of the current API and is important for debugging generated RHS code.

## Test Strategy

The migration should be driven by parity tests, not by ad hoc manual inspection.

### Test categories

#### 1. Backend contract tests

Purpose:

- ensure the public simulators behave the same regardless of backend implementation

Coverage:

- constructor behavior
- `set_problem_data`
- `set_tspan`
- `set_solver_parameters`
- `seed_rng`
- `transient`
- `shift_x0`
- `shift_tspan`

#### 2. Numerical parity tests

Purpose:

- compare C++ backend and PyOpenCL backend on the same models and parameters

Coverage:

- deterministic model, fixed stepper, single precision
- deterministic model, fixed stepper, double precision where supported
- adaptive stepper parity with tolerances
- stochastic model parity with fixed seed
- ensemble execution parity

Acceptance:

- exact match when that is realistic
- otherwise explicit tolerances by solver and precision

#### 3. Trajectory output tests

Purpose:

- verify time, state, derivative, auxiliary, and `n_stored` outputs match

Coverage:

- `nout` behavior
- `max_store` truncation behavior
- reshape semantics into `TrajectoryOutput`

#### 4. Feature output tests

Purpose:

- verify feature names and feature arrays match
- verify observer initialization and continuation behavior

Coverage:

- one-pass observers
- two-pass observers
- event timestamp storage
- observer reinitialization vs continuation

#### 5. Build-key and cache tests

Purpose:

- ensure code-affecting changes invalidate correctly

Coverage:

- changing RHS source changes the build key
- changing problem dimensions changes the build key
- changing observer or stored event count changes the build key
- repeated identical builds hit the cache

#### 6. Error-path tests

Purpose:

- ensure users get actionable failures

Coverage:

- malformed RHS
- unknown stepper
- unknown observer
- unsupported double precision on selected device

## Rollout Strategy

The backend should not switch from C++ to PyOpenCL in one change.

### Rollout mechanism

Add an internal backend selector, for example:

- default backend: `cpp`
- test backend: `pyopencl`

This selector can be driven by:

- an internal constructor argument not exposed in the public docs
- an environment variable for tests and local development
- a dedicated test helper

The public API should not advertise backend selection until the PyOpenCL path is ready.

### Exit condition for default switch

Do not make PyOpenCL the default until all of the following are true.

- transient parity is complete
- trajectory parity is complete
- feature parity is complete for the current shipped observer set
- error reporting is at least as good as the C++ path
- the packaging and dependency story is stable

## Phased Implementation Plan

The implementation should be shipped in small, reviewable phases.

### Phase 0: Backend seam and parity harness

Priority: P0

Goal:

- create the internal backend boundary without changing user-visible behavior

Deliverables:

- `clode/_backends/protocol.py`
- `clode/_backends/cpp.py`
- `clode/_backends/factory.py`
- public simulators instantiate backends through the factory
- parity test scaffolding can run the public API against the C++ backend through the new seam

Not in scope:

- PyOpenCL execution
- runtime rewrite

Acceptance criteria:

- existing public behavior unchanged
- no public API changes
- all existing tests pass

### Phase 1: PyOpenCL primitives

Priority: P0

Goal:

- build the reusable primitives without yet trying to support all solver modes

Deliverables:

- `clode/_pyopencl/models.py`
- `clode/_pyopencl/runtime.py`
- `clode/_pyopencl/registry.py`
- `clode/_pyopencl/source_builder.py`
- `clode/_pyopencl/program_cache.py`
- `clode/_pyopencl/errors.py`
- tests for build-key generation, source assembly, and build failures

Not in scope:

- trajectory execution
- feature execution

Acceptance criteria:

- runtime can compile a transient program
- `get_program_string()` equivalent source is available
- cache hit and miss behavior is tested

### Phase 2: Transient backend parity

Priority: P0

Goal:

- ship the first usable PyOpenCL execution path

Deliverables:

- `PyOpenCLTransientBackend`
- common buffer allocation and upload/download path
- transient kernel launch path
- parity tests against the C++ backend

Not in scope:

- trajectory buffers
- observers
- output model simplification

Acceptance criteria:

- deterministic transient parity on representative models
- stochastic parity with fixed seed where expected
- ensemble state transfer and continuation behavior match the C++ backend

### Phase 3: Trajectory backend parity

Priority: P1

Goal:

- add trajectory behavior without broad redesign

Deliverables:

- `PyOpenCLTrajectoryBackend`
- trajectory buffer allocation
- trajectory readback and reshape compatibility
- parity tests for `nout`, `max_store`, and `n_stored`

Not in scope:

- chunked streaming rewrite
- monitor abstraction

Acceptance criteria:

- `TrajectoryOutput` behavior unchanged
- trajectory arrays match current backend within agreed tolerances

### Phase 4: Feature backend parity

Priority: P1

Goal:

- support the current observer and feature pipeline as shipped today

Deliverables:

- `PyOpenCLFeatureBackend`
- observer parameter upload path
- observer warmup path
- feature buffer allocation
- observer continuation behavior
- parity tests for current shipped observers

Not in scope:

- observer redesign
- event storage redesign

Acceptance criteria:

- feature names match current backend
- feature arrays match current backend within agreed tolerances
- observer initialized vs reinitialized behavior matches current backend

### Phase 5: Controlled default switch

Priority: P1

Goal:

- make PyOpenCL the default backend only after parity is proven

Deliverables:

- backend selector defaults to PyOpenCL
- C++ backend remains available for fallback during transition
- docs updated for backend behavior and dependency requirements

Acceptance criteria:

- CI passes on supported platforms
- parity suite passes with PyOpenCL as default
- fallback path still works

### Phase 6: Post-parity simplification

Priority: P2

Goal:

- simplify architecture after the new backend is already shipping

Candidate tasks:

- unify duplicated kernel prolog and epilog logic
- revisit `transient`, `trajectory`, and `features` under a broader monitor or sink model
- separate persistent observer state from optional event-output storage
- add chunked trajectory streaming
- remove mixed `__cplusplus` metadata from OpenCL files
- retire the C++ backend when it is no longer needed

Acceptance criteria:

- no public contract regressions
- parity tests continue to pass
- complexity measurably decreases rather than merely moving around

## Risks And Mitigations

### Risk: source assembly diverges from current behavior

Mitigation:

- preserve current file composition in phase one
- expose `get_program_string()` for both backends
- compare generated source and build options in tests

### Risk: backend seam becomes too abstract

Mitigation:

- keep the protocol close to the current pybind shape
- avoid introducing generalized execution graphs or plugin systems early

### Risk: PyOpenCL dependency and packaging complexity

Mitigation:

- keep the C++ backend available during rollout
- add dependency and packaging work only after transient parity is working locally and in CI

### Risk: parity work gets derailed by early conceptual cleanup

Mitigation:

- explicitly defer monitor or sink unification to phase six
- review migration changes against preservation rules

## Deferred Design: Monitor Or Sink Abstraction

It is very likely that `transient`, `trajectory`, and `features` are special cases of a broader execution driver plus output policy model.

That future design probably looks like:

- base solver state owned by a shared driver
- optional warmup stage
- optional per-step sampling sink
- optional event sink
- optional reduction sink

However, that design should be deferred until after the PyOpenCL backend is already shipping and parity is established. It is the correct next simplification, not the correct first migration.

## Immediate Next Steps

1. Land the backend seam with a C++ adapter.
2. Add the PyOpenCL core models, source builder, and program cache.
3. Implement transient parity first.

Those three steps create a disciplined path to a shippable backend without forcing early redesign in the most coupled part of the current system.
