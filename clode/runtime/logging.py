from __future__ import annotations

from enum import IntEnum


class LogLevel(IntEnum):
    """Runtime logging levels used by clODE."""

    trace = 0
    debug = 1
    info = 2
    warn = 3
    err = 4
    critical = 5
    off = 6


DEFAULT_LOG_LEVEL = LogLevel.warn
_LOG_LEVEL = DEFAULT_LOG_LEVEL
_LOG_PATTERN: str | None = None


def _coerce_log_level(level: LogLevel | int) -> LogLevel:
    return LogLevel(int(level))


def get_log_level() -> LogLevel:
    return _LOG_LEVEL


def set_log_level(level: LogLevel | int) -> None:
    global _LOG_LEVEL

    _LOG_LEVEL = _coerce_log_level(level)


def set_log_pattern(pattern: str) -> None:
    global _LOG_PATTERN

    _LOG_PATTERN = pattern


__all__ = [
    "DEFAULT_LOG_LEVEL",
    "LogLevel",
    "get_log_level",
    "set_log_level",
    "set_log_pattern",
]