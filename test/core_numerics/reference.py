from __future__ import annotations

import math

import numpy as np


def _time_array(t: float | np.ndarray) -> np.ndarray:
    return np.asarray(t, dtype=np.float64)


def stable_linear_state(
    t: float | np.ndarray,
    *,
    x0: float = 2.0,
    y0: float = -1.5,
    a: float = 0.5,
    b: float = 1.25,
) -> np.ndarray:
    time = _time_array(t)
    x = x0 * np.exp(-a * time)
    y = y0 * np.exp(-b * time)
    return np.stack((x, y), axis=-1)


def stable_linear_derivative(
    t: float | np.ndarray,
    *,
    x0: float = 2.0,
    y0: float = -1.5,
    a: float = 0.5,
    b: float = 1.25,
) -> np.ndarray:
    state = stable_linear_state(t, x0=x0, y0=y0, a=a, b=b)
    x = state[..., 0]
    y = state[..., 1]
    return np.stack((-a * x, -b * y), axis=-1)


def stable_linear_aux(
    t: float | np.ndarray,
    *,
    x0: float = 2.0,
    y0: float = -1.5,
    a: float = 0.5,
    b: float = 1.25,
) -> np.ndarray:
    state = stable_linear_state(t, x0=x0, y0=y0, a=a, b=b)
    x = state[..., 0]
    y = state[..., 1]
    return np.stack((x + 2.0 * y, 3.0 * x - y), axis=-1)


def fixed_step_time_grid(start: float, end: float, dt: float) -> np.ndarray:
    current = np.float32(start)
    final = np.float32(end)
    step = np.float32(dt)
    times = [float(current)]

    while current < final:
        current = np.float32(current + step)
        times.append(float(current))

    return np.asarray(times, dtype=np.float64)


def fixed_step_observer_times(start: float, end: float, dt: float) -> np.ndarray:
    return fixed_step_time_grid(start, end, dt)[1:]


def fixed_step_step_count(start: float, end: float, dt: float) -> int:
    return int(fixed_step_observer_times(start, end, dt).size)


def fixed_step_stored_times(start: float, end: float, dt: float, nout: int) -> np.ndarray:
    return fixed_step_time_grid(start, end, dt)[::nout]


def hopf_on_cycle_state(t: float | np.ndarray, *, omega: float = 1.0) -> np.ndarray:
    time = _time_array(t)
    phase = omega * time
    return np.stack((np.cos(phase), np.sin(phase)), axis=-1)


def hopf_on_cycle_derivative(t: float | np.ndarray, *, omega: float = 1.0) -> np.ndarray:
    time = _time_array(t)
    phase = omega * time
    return np.stack((-omega * np.sin(phase), omega * np.cos(phase)), axis=-1)


def hopf_on_cycle_aux(t: float | np.ndarray) -> np.ndarray:
    time = _time_array(t)
    return np.ones(time.shape + (1,), dtype=np.float64)


def hopf_cycle_period(*, omega: float = 1.0) -> float:
    return 2.0 * math.pi / omega


def ou_stationary_mean(mu: float) -> float:
    return float(mu)


def ou_stationary_variance(sigma: float) -> float:
    return float(sigma * sigma / 2.0)
