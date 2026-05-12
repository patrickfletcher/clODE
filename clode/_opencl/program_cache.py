from __future__ import annotations

import logging
import os
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


LOGGER = logging.getLogger(__name__)


def _compiler_output_hint() -> str | None:
    if os.environ.get("PYOPENCL_COMPILER_OUTPUT"):
        return None
    return (
        "PyOpenCL compiler messages are hidden. Set PYOPENCL_COMPILER_OUTPUT=1 or "
        "call clode.configure_logging(..., compiler_output=True) before building to "
        "surface vendor compiler output."
    )


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
            LOGGER.debug("Reusing cached OpenCL program for build key %r", source_bundle.build_key)
            return existing

        opencl_binding = _require_opencl_binding()
        program = opencl_binding.Program(runtime.context, source_bundle.source_text)
        try:
            LOGGER.debug(
                "Building OpenCL program for build key %r with options=%s",
                source_bundle.build_key,
                list(source_bundle.build_options),
            )
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
        LOGGER.debug("Built OpenCL program for build key %r", source_bundle.build_key)
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
        compiler_output_hint = _compiler_output_hint()
        if compiler_output_hint is not None:
            details.append(compiler_output_hint)
        return "\n\n".join(detail for detail in details if detail)
