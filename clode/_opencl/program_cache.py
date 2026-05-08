from __future__ import annotations

from typing import TYPE_CHECKING, Any
import warnings

from .errors import BuildError, OpenCLDependencyError
from .models import ProgramBundle, SourceBundle

if TYPE_CHECKING:
    from .runtime import OpenCLRuntime

try:
    import pyopencl as cl
except ModuleNotFoundError:
    cl = None


def _require_opencl_binding() -> Any:
    if cl is None:
        raise OpenCLDependencyError(
            "pyopencl is required for clODE's OpenCL runtime. Install pyopencl into the active environment or reinstall clode with its default dependencies."
        )
    return cl


class ProgramCache:
    def __init__(self) -> None:
        self._cache: dict[object, ProgramBundle] = {}

    def __len__(self) -> int:
        return len(self._cache)

    def clear(self) -> None:
        self._cache.clear()

    def get_or_build(
        self, runtime: OpenCLRuntime, source_bundle: SourceBundle
    ) -> ProgramBundle:
        existing = self._cache.get(source_bundle.build_key)
        if existing is not None:
            return existing

        opencl_binding = _require_opencl_binding()
        program = opencl_binding.Program(runtime.context, source_bundle.source_text)
        try:
            program.build(options=list(source_bundle.build_options))
            kernels = {
                kernel_name: opencl_binding.Kernel(program, kernel_name)
                for kernel_name in source_bundle.kernel_names
            }
        except Exception as exc:
            raise BuildError(
                "OpenCL program build failed",
                source_text=source_bundle.source_text,
                build_options=source_bundle.build_options,
                build_log=self._extract_build_log(program, runtime, exc),
            ) from exc

        bundle = ProgramBundle(
            build_key=source_bundle.build_key,
            source_bundle=source_bundle,
            program=program,
            kernels=kernels,
        )
        self._cache[source_bundle.build_key] = bundle
        return bundle

    def _extract_build_log(
        self, program: object, runtime: OpenCLRuntime, exc: Exception
    ) -> str:
        opencl_binding = _require_opencl_binding()
        details: list[str] = []
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                build_log = program.get_build_info(
                    runtime.device, opencl_binding.program_build_info.LOG
                )
            if build_log:
                details.append(str(build_log))
        except Exception:
            pass
        if not details:
            details.append(str(exc))
        return "\n\n".join(detail for detail in details if detail)
