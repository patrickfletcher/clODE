from __future__ import annotations

import hashlib
from pathlib import Path
import re

from .errors import RhsValidationError, UnsupportedObserverError
from .models import (
    OPENCL_BACKEND_VERSION,
    BuildKey,
    KernelKind,
    Precision,
    ProblemShape,
    RhsSource,
    SourceBundle,
)
from .registry import KernelRegistry


class SourceBuilder:
    _GET_RHS_PATTERN = re.compile(r"\bvoid\s+getRHS\s*\(")

    def __init__(self, kernel_root: Path, registry: KernelRegistry | None = None) -> None:
        self._kernel_root = kernel_root
        self._registry = registry or KernelRegistry(kernel_root)

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
        self._validate_kernel_configuration(kernel_kind, observer_name, n_store_events)
        self._registry.validate_observer(observer_name)
        self._validate_rhs(rhs)

        entrypoint_paths = self._registry.get_entrypoint_paths(kernel_kind)
        entrypoint_source = "".join(
            path.read_text(encoding="utf-8") for path in entrypoint_paths
        )
        kernel_tree_digest = self._compute_kernel_tree_digest()

        build_key = BuildKey(
            backend_version=OPENCL_BACKEND_VERSION,
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

        return SourceBundle(
            build_key=build_key,
            source_text=entrypoint_source + rhs.text,
            build_options=self._make_build_options(build_key),
            kernel_names=self._registry.get_kernel_names(kernel_kind),
        )

    def _compute_kernel_tree_digest(self) -> str:
        digest = hashlib.sha256()
        for path in sorted(self._kernel_root.rglob("*.cl")) + sorted(
            self._kernel_root.rglob("*.clh")
        ):
            digest.update(path.relative_to(self._kernel_root).as_posix().encode("utf-8"))
            digest.update(b"\0")
            digest.update(path.read_bytes())
            digest.update(b"\0")
        return digest.hexdigest()

    def _make_build_options(self, key: BuildKey) -> tuple[str, ...]:
        options = [
            "-DCLODE_SINGLE_PRECISION"
            if key.precision is Precision.SINGLE
            else "-DCLODE_DOUBLE_PRECISION",
            f"-D{self._registry.get_stepper_define(key.stepper_name)}",
            f"-DN_PAR={key.problem_shape.n_par}",
            f"-DN_VAR={key.problem_shape.n_var}",
            f"-DN_AUX={key.problem_shape.n_aux}",
            f"-DN_WIENER={key.problem_shape.n_wiener}",
            f"-I{self._kernel_root}",
        ]
        if key.observer_name is not None:
            options.append(f"-D{self._registry.get_observer_define(key.observer_name)}")
            options.append(f"-DN_STORE_EVENTS={key.n_store_events}")
        if key.debug_build:
            options.extend(["-g", "-cl-opt-disable"])
        return tuple(options)

    def _validate_kernel_configuration(
        self,
        kernel_kind: KernelKind,
        observer_name: str | None,
        n_store_events: int,
    ) -> None:
        if kernel_kind is KernelKind.FEATURES:
            if observer_name is None:
                raise UnsupportedObserverError(observer_name)
            return
        if observer_name is not None:
            raise ValueError("observer_name is only valid for features builds")
        if n_store_events != 0:
            raise ValueError("n_store_events is only valid for features builds")

    def _validate_rhs(self, rhs: RhsSource) -> None:
        if not self._GET_RHS_PATTERN.search(rhs.text):
            raise RhsValidationError(
                "RHS source must define void getRHS(...)",
                origin_label=rhs.origin_label,
            )
