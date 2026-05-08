from __future__ import annotations


class OpenCLBackendError(RuntimeError):
    pass


class OpenCLDependencyError(OpenCLBackendError):
    pass


class OpenCLValidationError(OpenCLBackendError, ValueError):
    pass


class RhsValidationError(OpenCLValidationError):
    def __init__(self, message: str, origin_label: str | None = None) -> None:
        detail = message if origin_label is None else f"{message} (rhs={origin_label})"
        super().__init__(detail)
        self.origin_label = origin_label


class RegistryValidationError(OpenCLValidationError):
    pass


class UnsupportedStepperError(RegistryValidationError):
    def __init__(self, stepper_name: str) -> None:
        super().__init__(f"Unsupported stepper name: {stepper_name!r}")
        self.stepper_name = stepper_name


class UnsupportedObserverError(RegistryValidationError):
    def __init__(self, observer_name: str | None) -> None:
        super().__init__(f"Unsupported observer name: {observer_name!r}")
        self.observer_name = observer_name


class DoublePrecisionNotSupportedError(OpenCLBackendError):
    def __init__(self, device_name: str | None = None) -> None:
        detail = "Selected device does not support double precision"
        if device_name is not None:
            detail = f"{detail}: {device_name}"
        super().__init__(detail)
        self.device_name = device_name


class BuildError(OpenCLBackendError):
    def __init__(
        self,
        message: str,
        *,
        source_text: str,
        build_options: tuple[str, ...],
        build_log: str,
    ) -> None:
        super().__init__(message)
        self.source_text = source_text
        self.build_options = build_options
        self.build_log = build_log

    @property
    def formatted_build_options(self) -> str:
        return " ".join(self.build_options)
