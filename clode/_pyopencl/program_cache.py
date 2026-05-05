from __future__ import annotations

from typing import TYPE_CHECKING, Any
import warnings

from .errors import BuildError, PyOpenCLDependencyError
from .models import ProgramBundle, SourceBundle

if TYPE_CHECKING:
    from .runtime import OpenCLRuntime

try:
    import pyopencl as cl
except ModuleNotFoundError:
    cl = None


def _require_pyopencl() -> Any:
    if cl is None:
        raise PyOpenCLDependencyError(
            "pyopencl is required for the PyOpenCL backend. Install pyopencl directly or use the optional 'clode[pyopencl]' dependency."
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

        pyopencl = _require_pyopencl()
        program = pyopencl.Program(runtime.context, source_bundle.source_text)
        try:
            program.build(options=list(source_bundle.build_options))
            kernels = {
                kernel_name: pyopencl.Kernel(program, kernel_name)
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
        pyopencl = _require_pyopencl()
        details: list[str] = []
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", UserWarning)
                build_log = program.get_build_info(
                    runtime.device, pyopencl.program_build_info.LOG
                )
            if build_log:
                details.append(str(build_log))
        except Exception:
            pass
        if not details:
            details.append(str(exc))
        return "\n\n".join(detail for detail in details if detail)
