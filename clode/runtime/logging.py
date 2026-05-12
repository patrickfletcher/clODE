from __future__ import annotations

from collections.abc import Sequence
import logging as pylogging
import os
from typing import TextIO


_LOGGER_NAME = "clode"
_PYWARNINGS_LOGGER_NAME = "py.warnings"
_DEFAULT_LOG_LEVEL = pylogging.WARNING
_DEFAULT_LOG_FORMAT = "%(name)s %(levelname)s: %(message)s"


def _normalize_level(level: int | str) -> int:
    if isinstance(level, str):
        normalized = pylogging.getLevelNamesMapping().get(level.upper())
        if normalized is None:
            raise ValueError(f"Unknown logging level: {level!r}")
        return int(normalized)
    return int(level)


def _replace_handlers(logger: pylogging.Logger, handlers: Sequence[pylogging.Handler]) -> None:
    for existing in list(logger.handlers):
        logger.removeHandler(existing)
        existing.close()
    for handler in handlers:
        logger.addHandler(handler)


def _make_stream_handler(
    stream: TextIO | None,
    level: int,
    fmt: str,
) -> pylogging.Handler:
    handler = pylogging.StreamHandler(stream)
    handler.setLevel(level)
    handler.setFormatter(pylogging.Formatter(fmt))
    return handler


def _set_compiler_output(enabled: bool | None) -> None:
    if enabled is None:
        return
    if enabled:
        os.environ["PYOPENCL_COMPILER_OUTPUT"] = "1"
    else:
        os.environ.pop("PYOPENCL_COMPILER_OUTPUT", None)


def get_logger(name: str | None = None) -> pylogging.Logger:
    """Return the root clODE logger or one of its child loggers."""

    if not name:
        logger_name = _LOGGER_NAME
    elif name == _LOGGER_NAME or name.startswith(f"{_LOGGER_NAME}."):
        logger_name = name
    else:
        logger_name = f"{_LOGGER_NAME}.{name}"

    logger = pylogging.getLogger(logger_name)
    if logger_name == _LOGGER_NAME and not logger.handlers:
        logger.addHandler(pylogging.NullHandler())
    return logger


def configure_logging(
    level: int | str = _DEFAULT_LOG_LEVEL,
    *,
    stream: TextIO | None = None,
    format: str = _DEFAULT_LOG_FORMAT,
    capture_warnings: bool = True,
    compiler_output: bool | None = None,
) -> pylogging.Logger:
    """Configure clODE logging and optional PyOpenCL diagnostics."""

    normalized_level = _normalize_level(level)
    logger = get_logger()
    handler = _make_stream_handler(stream, normalized_level, format)
    _replace_handlers(logger, [handler])
    logger.setLevel(normalized_level)
    logger.propagate = False

    warning_logger = pylogging.getLogger(_PYWARNINGS_LOGGER_NAME)
    if capture_warnings:
        pylogging.captureWarnings(True)
        warning_handler = _make_stream_handler(stream, pylogging.WARNING, format)
        _replace_handlers(warning_logger, [warning_handler])
        warning_logger.setLevel(pylogging.WARNING)
        warning_logger.propagate = False
    else:
        pylogging.captureWarnings(False)
        _replace_handlers(warning_logger, [])

    _set_compiler_output(compiler_output)
    logger.debug(
        "Configured clODE logging level=%s capture_warnings=%s compiler_output=%s",
        pylogging.getLevelName(normalized_level),
        capture_warnings,
        compiler_output,
    )
    return logger


get_logger().setLevel(_DEFAULT_LOG_LEVEL)


__all__ = [
    "configure_logging",
    "get_logger",
]