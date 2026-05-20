from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# This public example mirrors the kernel formulas in NumPy so the tradeoffs are
# easy to inspect and rerun outside OpenCL debugging. The actual OpenCL helper
# prototypes are pinned separately in test/kernel_components/test_kernel_math.py.

try:
    import matplotlib.pyplot as plt
except ImportError:  # pragma: no cover - plotting is optional for the example
    plt = None


@dataclass(frozen=True)
class MeanScenarioResult:
    name: str
    expected: float
    running_mean_time: float
    compensated_mean: float

    @property
    def running_error(self) -> float:
        return self.running_mean_time - self.expected

    @property
    def compensated_error(self) -> float:
        return self.compensated_mean - self.expected


@dataclass(frozen=True)
class TimeMethodResult:
    method: str
    final_time: float
    final_error: float
    max_error: float
    first_stalled_step: int | None
    unique_increment_count: int


@dataclass(frozen=True)
class TimeScenarioResult:
    name: str
    t0: float
    dt: float
    n_steps: int
    spacing_at_t0: float
    methods: tuple[TimeMethodResult, ...]


@dataclass(frozen=True)
class MeanOriginScenarioResult:
    t0: float
    expected: float
    elapsed_difference_mean: float
    relative_elapsed_mean: float
    elapsed_difference_elapsed: float
    relative_elapsed: float


@dataclass(frozen=True)
class ObserverGuardScenarioResult:
    name: str
    baseline_label: str
    baseline_event_count: int
    guard_label: str
    guarded_event_count: int
    note: str


@dataclass(frozen=True)
class LocalExtremumPrototypeResult:
    name: str
    argmax_mean_time_error: float
    argmax_max_time_error: float
    quadratic_mean_time_error: float
    quadratic_max_time_error: float
    argmax_mean_value_error: float
    quadratic_mean_value_error: float
    note: str


@dataclass(frozen=True)
class ThresholdInterpolationPrototypeResult:
    name: str
    sample_mean_time_error: float
    sample_max_time_error: float
    linear_mean_time_error: float
    linear_max_time_error: float
    hermite_mean_time_error: float
    hermite_max_time_error: float
    note: str


STEP_COUNTER_CHUNK_SIZE = 1_000_000
LONG_HORIZON_STEPS = 30_000_000


def running_mean_time_float32(values: np.ndarray, dt: float = 1.0) -> np.float32:
    mean = np.float32(0.0)
    total_delta = np.float32(0.0)
    dt32 = np.float32(dt)
    for value in values.astype(np.float32, copy=False):
        total_delta = np.float32(total_delta + dt32)
        update = np.float32((value - mean) * np.float32(dt32 / total_delta))
        mean = np.float32(mean + update)
    return mean


def compensated_integral_mean_float32(values: np.ndarray, dt: float = 1.0) -> np.float32:
    integral = np.float32(0.0)
    correction = np.float32(0.0)
    total_delta = np.float32(0.0)
    dt32 = np.float32(dt)
    for value in values.astype(np.float32, copy=False):
        total_delta = np.float32(total_delta + dt32)
        new_value = np.float32(dt32 * value)
        total = np.float32(integral + new_value)
        if abs(float(integral)) >= abs(float(new_value)):
            correction = np.float32(correction + np.float32((integral - total) + new_value))
        else:
            correction = np.float32(correction + np.float32((new_value - total) + integral))
        integral = total
    return np.float32((integral + correction) / total_delta)


def compensated_mean_with_elapsed_difference_float32(
    values: np.ndarray,
    *,
    t0: float,
    dt: float,
) -> tuple[float, float]:
    integral = np.float32(0.0)
    correction = np.float32(0.0)
    t_start = np.float32(t0)
    ti = np.float32(t0)
    dt32 = np.float32(dt)
    mean_value = 0.0
    elapsed = np.float32(0.0)
    for value in values.astype(np.float32, copy=False):
        ti = np.float32(ti + dt32)
        new_value = np.float32(dt32 * value)
        total = np.float32(integral + new_value)
        if abs(float(integral)) >= abs(float(new_value)):
            correction = np.float32(correction + np.float32((integral - total) + new_value))
        else:
            correction = np.float32(correction + np.float32((new_value - total) + integral))
        integral = total
        elapsed = np.float32(ti - t_start)
        if elapsed > 0.0:
            mean_value = float(np.float32((integral + correction) / elapsed))
    return mean_value, float(elapsed)


def compensated_mean_with_relative_elapsed_float32(
    values: np.ndarray,
    *,
    dt: float,
) -> tuple[float, float]:
    integral = np.float32(0.0)
    correction = np.float32(0.0)
    elapsed = np.float32(0.0)
    dt32 = np.float32(dt)
    mean_value = 0.0
    for value in values.astype(np.float32, copy=False):
        elapsed = np.float32(elapsed + dt32)
        new_value = np.float32(dt32 * value)
        total = np.float32(integral + new_value)
        if abs(float(integral)) >= abs(float(new_value)):
            correction = np.float32(correction + np.float32((integral - total) + new_value))
        else:
            correction = np.float32(correction + np.float32((new_value - total) + integral))
        integral = total
        if elapsed > 0.0:
            mean_value = float(np.float32((integral + correction) / elapsed))
    return mean_value, float(elapsed)


def direct_time_updates_float32(t0: float, dt: float, n_steps: int) -> np.ndarray:
    times = np.empty(n_steps + 1, dtype=np.float32)
    dt32 = np.float32(dt)
    times[0] = np.float32(t0)
    for step in range(1, n_steps + 1):
        times[step] = np.float32(times[step - 1] + dt32)
    return times


def kahan_time_updates_float32(t0: float, dt: float, n_steps: int) -> np.ndarray:
    times = np.empty(n_steps + 1, dtype=np.float32)
    dt32 = np.float32(dt)
    correction = np.float32(0.0)
    times[0] = np.float32(t0)
    for step in range(1, n_steps + 1):
        y = np.float32(dt32 - correction)
        total = np.float32(times[step - 1] + y)
        correction = np.float32((total - times[step - 1]) - y)
        times[step] = total
    return times


def structured_time_updates_float32(t0: float, dt: float, n_steps: int) -> np.ndarray:
    steps = np.arange(n_steps + 1, dtype=np.float32)
    return np.float32(np.float32(t0) + steps * np.float32(dt))


def summarize_structured_step_counter_product_float32(
    method: str,
    t0: float,
    dt: float,
    n_steps: int,
    *,
    chunk_size: int = STEP_COUNTER_CHUNK_SIZE,
) -> TimeMethodResult:
    dt32 = np.float32(dt)
    t032 = np.float32(t0)
    max_error = 0.0
    previous_time: np.float32 | None = None
    first_stalled_step: int | None = None
    unique_increments: set[float] = set()
    final_time = float(t0)

    for start in range(0, n_steps + 1, chunk_size):
        stop = min(start + chunk_size, n_steps + 1)
        steps = np.arange(start, stop, dtype=np.int64)
        times = np.float32(t032 + steps.astype(np.float32) * dt32)
        reference = t0 + steps.astype(np.float64) * dt
        max_error = max(max_error, float(np.abs(times.astype(np.float64) - reference).max()))

        if previous_time is not None:
            boundary_increment = float(np.float32(times[0] - previous_time))
            unique_increments.add(boundary_increment)
            if first_stalled_step is None and boundary_increment == 0.0:
                first_stalled_step = start

        increments = np.diff(times)
        if increments.size:
            unique_increments.update(float(value) for value in np.unique(increments))
            if first_stalled_step is None:
                zero_steps = np.nonzero(increments == 0.0)[0]
                if zero_steps.size != 0:
                    first_stalled_step = start + int(zero_steps[0]) + 1

        previous_time = np.float32(times[-1])
        final_time = float(times[-1])

    return TimeMethodResult(
        method=method,
        final_time=final_time,
        final_error=abs(final_time - (t0 + n_steps * dt)),
        max_error=max_error,
        first_stalled_step=first_stalled_step,
        unique_increment_count=len(unique_increments),
    )


def summarize_direct_add_chunked_float32(
    method: str,
    t0: float,
    dt: float,
    n_steps: int,
    *,
    chunk_size: int = STEP_COUNTER_CHUNK_SIZE,
) -> TimeMethodResult:
    t = np.float32(t0)
    dt32 = np.float32(dt)
    first_stalled_step = None
    max_error = 0.0
    unique_increments: set[float] = set()

    for start in range(0, n_steps, chunk_size):
        size = min(chunk_size, n_steps - start)
        trace = np.empty(size + 1, dtype=np.float32)
        trace[0] = t
        trace[1:] = dt32
        trace = np.add.accumulate(trace, dtype=np.float32)
        increments = np.diff(trace)
        unique_increments.update(float(value) for value in np.unique(increments))
        if first_stalled_step is None:
            zero_steps = np.nonzero(increments == 0.0)[0]
            if zero_steps.size != 0:
                first_stalled_step = start + int(zero_steps[0]) + 1

        reference = t0 + np.arange(start + 1, start + size + 1, dtype=np.float64) * dt
        max_error = max(max_error, float(np.abs(trace[1:].astype(np.float64) - reference).max()))
        t = np.float32(trace[-1])

    final_time = float(t)
    return TimeMethodResult(
        method=method,
        final_time=final_time,
        final_error=abs(final_time - (t0 + n_steps * dt)),
        max_error=max_error,
        first_stalled_step=first_stalled_step,
        unique_increment_count=len(unique_increments),
    )


def summarize_kahan_add_float32(
    method: str,
    t0: float,
    dt: float,
    n_steps: int,
) -> TimeMethodResult:
    t = np.float32(t0)
    correction = np.float32(0.0)
    dt32 = np.float32(dt)
    previous = t
    first_stalled_step = None
    max_error = 0.0
    unique_increments: set[float] = set()

    for step in range(1, n_steps + 1):
        y = np.float32(dt32 - correction)
        total = np.float32(t + y)
        correction = np.float32((total - t) - y)
        t = total
        increment = float(np.float32(t - previous))
        unique_increments.add(increment)
        if first_stalled_step is None and increment == 0.0:
            first_stalled_step = step
        previous = t
        reference = t0 + step * dt
        max_error = max(max_error, abs(float(t) - reference))

    final_time = float(t)
    return TimeMethodResult(
        method=method,
        final_time=final_time,
        final_error=abs(final_time - (t0 + n_steps * dt)),
        max_error=max_error,
        first_stalled_step=first_stalled_step,
        unique_increment_count=len(unique_increments),
    )


def first_stalled_step(times: np.ndarray) -> int | None:
    zero_steps = np.nonzero(np.diff(times) == 0.0)[0]
    if zero_steps.size == 0:
        return None
    return int(zero_steps[0] + 1)


def summarize_time_method(
    method: str,
    times: np.ndarray,
    reference: np.ndarray,
) -> TimeMethodResult:
    error = np.abs(times.astype(np.float64) - reference)
    increments = np.diff(times)
    return TimeMethodResult(
        method=method,
        final_time=float(times[-1]),
        final_error=float(error[-1]),
        max_error=float(error.max()),
        first_stalled_step=first_stalled_step(times),
        unique_increment_count=int(np.unique(increments).size),
    )


def threshold2_style_event_count_float32(
    x: np.ndarray,
    dx: np.ndarray,
    *,
    min_amp: float = 0.0,
    x_up_threshold: float = 0.3,
    x_down_threshold: float = 0.2,
    dx_up_threshold: float = 0.0,
    dx_down_threshold: float = 0.0,
) -> int:
    x32 = np.asarray(x, dtype=np.float32)
    dx32 = np.asarray(dx, dtype=np.float32)

    x_global_max = float(np.max(x32))
    x_global_min = float(np.min(x32))
    dx_global_max = float(np.max(dx32))
    dx_global_min = float(np.min(dx32))

    amplitude = x_global_max - x_global_min
    x_up = x_global_min + x_up_threshold * amplitude
    x_down = x_global_min + (x_down_threshold * amplitude if x_down_threshold > 0.0 else x_up)
    dx_up = dx_up_threshold * dx_global_max
    dx_down = dx_down_threshold * dx_global_min if dx_down_threshold > 0.0 else dx_global_min

    if amplitude < min_amp:
        return 0

    in_upstate = bool(x32[0] > x_up)
    event_count = 0
    for xi, dxi in zip(x32[1:], dx32[1:]):
        xi_float = float(xi)
        dxi_float = float(dxi)
        if in_upstate:
            if xi_float <= x_down and dxi_float >= dx_down:
                in_upstate = False
        else:
            if xi_float > x_up and dxi_float > dx_up:
                event_count += 1
                in_upstate = True
    return event_count


def linear_crossing_time(
    t0: float,
    t1: float,
    x0: float,
    x1: float,
    x_target: float,
) -> float:
    return t0 + (x_target - x0) * (t1 - t0) / (x1 - x0)


def hermite_value_and_derivative(
    u: float,
    x0: float,
    x1: float,
    dx0: float,
    dx1: float,
    dt: float,
) -> tuple[float, float]:
    h00 = 2.0 * u**3 - 3.0 * u**2 + 1.0
    h10 = u**3 - 2.0 * u**2 + u
    h01 = -2.0 * u**3 + 3.0 * u**2
    h11 = u**3 - u**2
    value = h00 * x0 + h10 * dt * dx0 + h01 * x1 + h11 * dt * dx1

    dh00 = 6.0 * u**2 - 6.0 * u
    dh10 = 3.0 * u**2 - 4.0 * u + 1.0
    dh01 = -6.0 * u**2 + 6.0 * u
    dh11 = 3.0 * u**2 - 2.0 * u
    derivative = dh00 * x0 + dh10 * dt * dx0 + dh01 * x1 + dh11 * dt * dx1
    return value, derivative


def hermite_crossing_time(
    t0: float,
    t1: float,
    x0: float,
    x1: float,
    dx0: float,
    dx1: float,
    x_target: float,
    *,
    iterations: int = 8,
) -> float:
    dt = t1 - t0
    u = np.clip((x_target - x0) / (x1 - x0), 0.0, 1.0)
    for _ in range(iterations):
        value, derivative = hermite_value_and_derivative(u, x0, x1, dx0, dx1, dt)
        if abs(derivative) < 1.0e-12:
            break
        candidate = u - (value - x_target) / derivative
        if candidate <= 0.0 or candidate >= 1.0:
            break
        u = candidate
    return t0 + u * dt


def build_mean_scenarios() -> tuple[MeanScenarioResult, ...]:
    step_offset_values = np.ones(200_000, dtype=np.float64)
    step_offset_values[100_000:] += 1.0e-4

    oscillation_parameter = np.linspace(0.0, 1_000.0, 200_000, dtype=np.float64)
    smooth_oscillation_values = 1.0 + 1.0e-4 * np.sin(oscillation_parameter)

    scenarios = (
        ("half-window 1e-4 offset", step_offset_values),
        ("smooth 1 +/- 1e-4 sin(t)", smooth_oscillation_values),
    )

    results: list[MeanScenarioResult] = []
    for name, values in scenarios:
        results.append(
            MeanScenarioResult(
                name=name,
                expected=float(np.mean(values, dtype=np.float64)),
                running_mean_time=float(running_mean_time_float32(values)),
                compensated_mean=float(compensated_integral_mean_float32(values)),
            )
        )
    return tuple(results)


def build_mean_origin_scenarios() -> tuple[MeanOriginScenarioResult, ...]:
    values = np.ones(10_000, dtype=np.float64)
    values[5_000:] += 1.0e-4
    expected = float(np.mean(values, dtype=np.float64))
    relative_elapsed_mean, relative_elapsed = compensated_mean_with_relative_elapsed_float32(
        values,
        dt=0.01,
    )

    results: list[MeanOriginScenarioResult] = []
    for t0 in (0.0, 1.0e4, 1.0e6):
        elapsed_difference_mean, elapsed_difference_elapsed = (
            compensated_mean_with_elapsed_difference_float32(
                values,
                t0=t0,
                dt=0.01,
            )
        )
        results.append(
            MeanOriginScenarioResult(
                t0=t0,
                expected=expected,
                elapsed_difference_mean=elapsed_difference_mean,
                relative_elapsed_mean=relative_elapsed_mean,
                elapsed_difference_elapsed=elapsed_difference_elapsed,
                relative_elapsed=relative_elapsed,
            )
        )
    return tuple(results)


def build_time_scenarios() -> tuple[TimeScenarioResult, ...]:
    scenario_definitions = (
        ("zero-origin drift", 0.0, 0.01, 200_000),
        ("large-origin stall", 1.0e6, 0.01, 10_000),
    )

    results: list[TimeScenarioResult] = []
    for name, t0, dt, n_steps in scenario_definitions:
        reference = t0 + np.arange(n_steps + 1, dtype=np.float64) * dt
        methods = (
            summarize_time_method("direct add", direct_time_updates_float32(t0, dt, n_steps), reference),
            summarize_time_method("Kahan add", kahan_time_updates_float32(t0, dt, n_steps), reference),
            summarize_time_method(
                "structured t0 + step*dt",
                structured_time_updates_float32(t0, dt, n_steps),
                reference,
            ),
        )
        results.append(
            TimeScenarioResult(
                name=name,
                t0=t0,
                dt=dt,
                n_steps=n_steps,
                spacing_at_t0=float(np.spacing(np.float32(t0))),
                methods=methods,
            )
        )

    results.append(
        TimeScenarioResult(
            name="long zero-origin horizon",
            t0=0.0,
            dt=0.01,
            n_steps=LONG_HORIZON_STEPS,
            spacing_at_t0=float(np.spacing(np.float32(0.0))),
            methods=(
                summarize_direct_add_chunked_float32(
                    "direct add",
                    0.0,
                    0.01,
                    LONG_HORIZON_STEPS,
                ),
                summarize_kahan_add_float32(
                    "Kahan add",
                    0.0,
                    0.01,
                    LONG_HORIZON_STEPS,
                ),
                summarize_structured_step_counter_product_float32(
                    "structured t0 + step*dt",
                    0.0,
                    0.01,
                    LONG_HORIZON_STEPS,
                ),
            ),
        )
    )
    return tuple(results)


def build_observer_guard_scenarios() -> tuple[ObserverGuardScenarioResult, ...]:
    omega = 1.0
    sample_times = np.linspace(0.0, 40.0 * np.pi, 40_000, dtype=np.float64)

    tiny_amplitude = 5.0e-4
    tiny_x = tiny_amplitude * np.sin(omega * sample_times)
    tiny_dx = tiny_amplitude * omega * np.cos(omega * sample_times)

    ripple_x = np.sin(sample_times) + 0.1 * np.sin(55.0 * sample_times)
    ripple_dx = np.cos(sample_times) + 0.1 * 55.0 * np.cos(55.0 * sample_times)

    slope_gate_x = np.sin(sample_times) + 0.85 * np.sin(5.0 * sample_times)
    slope_gate_dx = np.cos(sample_times) + 0.85 * 5.0 * np.cos(5.0 * sample_times)

    return (
        ObserverGuardScenarioResult(
            name="tiny oscillation amplitude gate",
            baseline_label="min_amp = 0",
            baseline_event_count=threshold2_style_event_count_float32(
                tiny_x,
                tiny_dx,
                min_amp=0.0,
            ),
            guard_label="min_amp = 2e-3",
            guarded_event_count=threshold2_style_event_count_float32(
                tiny_x,
                tiny_dx,
                min_amp=2.0e-3,
            ),
            note="A trajectory range of about 1e-3 is treated as oscillatory without a floor, but suppressed once the minimum meaningful amplitude is set above that range.",
        ),
        ObserverGuardScenarioResult(
            name="threshold chatter with ripple",
            baseline_label="x_up = x_down = 0.5",
            baseline_event_count=threshold2_style_event_count_float32(
                ripple_x,
                ripple_dx,
                x_up_threshold=0.5,
                x_down_threshold=0.5,
            ),
            guard_label="x_up = 0.65, x_down = 0.35",
            guarded_event_count=threshold2_style_event_count_float32(
                ripple_x,
                ripple_dx,
                x_up_threshold=0.65,
                x_down_threshold=0.35,
            ),
            note="A single threshold counts many ripple-driven crossings, while a clear Schmitt-trigger-style hysteresis gap returns the expected coarse-cycle count.",
        ),
        ObserverGuardScenarioResult(
            name="slope gate on noisy crossings",
            baseline_label="dx_up = 0",
            baseline_event_count=threshold2_style_event_count_float32(
                slope_gate_x,
                slope_gate_dx,
                x_up_threshold=0.5,
                x_down_threshold=0.5,
                dx_up_threshold=0.0,
            ),
            guard_label="dx_up = 0.9",
            guarded_event_count=threshold2_style_event_count_float32(
                slope_gate_x,
                slope_gate_dx,
                x_up_threshold=0.5,
                x_down_threshold=0.5,
                dx_up_threshold=0.9,
            ),
            note="A derivative gate can suppress noisy shallow crossings when value thresholds alone are not selective enough.",
        ),
    )


def quadratic_vertex_from_three_samples(t: np.ndarray, y: np.ndarray) -> tuple[float, float]:
    b0 = y[0]
    b1 = (y[1] - b0) / (t[1] - t[0])
    b2 = (y[2] - b0 - b1 * (t[2] - t[0])) / ((t[2] - t[0]) * (t[2] - t[1]))
    tv = -(b1 - b2 * (t[0] + t[1])) / (2.0 * b2)
    yv = b0 + b1 * (tv - t[0]) + b2 * (tv - t[0]) * (tv - t[1])
    return float(tv), float(yv)


def build_local_extremum_prototype() -> LocalExtremumPrototypeResult:
    dt = 0.2
    offset = 0.03
    sample_times = np.arange(0.0, 20.0 * np.pi, dt, dtype=np.float64) + offset
    x = np.sin(sample_times)
    dx = np.cos(sample_times)

    argmax_time_errors: list[float] = []
    argmax_value_errors: list[float] = []
    quadratic_time_errors: list[float] = []
    quadratic_value_errors: list[float] = []

    for idx in range(2, len(sample_times)):
        if not (dx[idx - 1] > 0.0 and dx[idx] < 0.0):
            continue
        time_buffer = sample_times[idx - 2 : idx + 1]
        value_buffer = x[idx - 2 : idx + 1]
        true_peak_time = np.pi / 2.0 + 2.0 * np.pi * np.round(
            (time_buffer[1] - np.pi / 2.0) / (2.0 * np.pi)
        )
        if not (time_buffer[0] <= true_peak_time <= time_buffer[2]):
            continue

        argmax_index = int(np.argmax(value_buffer))
        argmax_time = float(time_buffer[argmax_index])
        argmax_value = float(value_buffer[argmax_index])
        quadratic_time, quadratic_value = quadratic_vertex_from_three_samples(
            time_buffer,
            value_buffer,
        )

        argmax_time_errors.append(abs(argmax_time - true_peak_time))
        argmax_value_errors.append(abs(argmax_value - 1.0))
        quadratic_time_errors.append(abs(quadratic_time - true_peak_time))
        quadratic_value_errors.append(abs(quadratic_value - 1.0))

    return LocalExtremumPrototypeResult(
        name="three-sample local maximum",
        argmax_mean_time_error=float(np.mean(argmax_time_errors)),
        argmax_max_time_error=float(np.max(argmax_time_errors)),
        quadratic_mean_time_error=float(np.mean(quadratic_time_errors)),
        quadratic_max_time_error=float(np.max(quadratic_time_errors)),
        argmax_mean_value_error=float(np.mean(argmax_value_errors)),
        quadratic_mean_value_error=float(np.mean(quadratic_value_errors)),
        note="On a coarse sampled sine wave, a quadratic vertex estimate on the same three-sample buffer is much more accurate than taking the buffered sample maximum directly.",
    )


def build_threshold_interpolation_prototype() -> ThresholdInterpolationPrototypeResult:
    dt = 0.2
    offset = 0.03
    threshold = 0.5
    sample_times = np.arange(0.0, 20.0 * 2.0 * np.pi, dt, dtype=np.float64) + offset
    x = np.sin(sample_times)
    dx = np.cos(sample_times)

    sample_errors: list[float] = []
    linear_errors: list[float] = []
    hermite_errors: list[float] = []

    base_crossing = float(np.arcsin(threshold))
    for idx in range(1, len(sample_times)):
        if not (x[idx - 1] <= threshold < x[idx]):
            continue
        true_time = base_crossing + 2.0 * np.pi * np.round(
            (sample_times[idx - 1] - base_crossing) / (2.0 * np.pi)
        )
        if true_time < sample_times[idx - 1]:
            true_time += 2.0 * np.pi
        if not (sample_times[idx - 1] <= true_time <= sample_times[idx]):
            continue

        sample_errors.append(abs(sample_times[idx] - true_time))
        linear_errors.append(
            abs(
                linear_crossing_time(
                    sample_times[idx - 1],
                    sample_times[idx],
                    x[idx - 1],
                    x[idx],
                    threshold,
                )
                - true_time
            )
        )
        hermite_errors.append(
            abs(
                hermite_crossing_time(
                    sample_times[idx - 1],
                    sample_times[idx],
                    x[idx - 1],
                    x[idx],
                    dx[idx - 1],
                    dx[idx],
                    threshold,
                )
                - true_time
            )
        )

    return ThresholdInterpolationPrototypeResult(
        name="upward threshold crossing on coarse sine",
        sample_mean_time_error=float(np.mean(sample_errors)),
        sample_max_time_error=float(np.max(sample_errors)),
        linear_mean_time_error=float(np.mean(linear_errors)),
        linear_max_time_error=float(np.max(linear_errors)),
        hermite_mean_time_error=float(np.mean(hermite_errors)),
        hermite_max_time_error=float(np.max(hermite_errors)),
        note="Linear inversion already improves coarse threshold timestamps substantially; slope-aware Hermite interpolation is even more accurate on smooth monotone crossings but needs stronger robustness checks before broad kernel use.",
    )


def print_mean_summary(results: tuple[MeanScenarioResult, ...]) -> None:
    print("Single-precision time-weighted means")
    print(
        "scenario                              exact            runningMeanTime   compensatedIntegral   |err running|   |err comp|"
    )
    for result in results:
        print(
            f"{result.name:<35}"
            f" {result.expected:>13.10f}"
            f" {result.running_mean_time:>18.10f}"
            f" {result.compensated_mean:>21.10f}"
            f" {abs(result.running_error):>15.3e}"
            f" {abs(result.compensated_error):>12.3e}"
        )
    print()


def print_mean_origin_summary(results: tuple[MeanOriginScenarioResult, ...]) -> None:
    print("Feature-window origin and elapsed-time denominator")
    print(
        "t0               exact            mean via ti-t_start   mean via relative elapsed   elapsed(ti-t_start)"
    )
    for result in results:
        print(
            f"{result.t0:>10.1f}"
            f" {result.expected:>13.10f}"
            f" {result.elapsed_difference_mean:>22.10f}"
            f" {result.relative_elapsed_mean:>26.10f}"
            f" {result.elapsed_difference_elapsed:>21.10f}"
        )
    print()


def print_time_summary(results: tuple[TimeScenarioResult, ...]) -> None:
    print("Single-precision time accumulation")
    for scenario in results:
        print(
            f"{scenario.name}: t0={scenario.t0:g}, dt={scenario.dt:g}, steps={scenario.n_steps}, "
            f"ulp(t0)={scenario.spacing_at_t0:g}"
        )
        print(
            "method                    final time       final err       max err    first stalled step   unique increments"
        )
        for method in scenario.methods:
            stalled = "none" if method.first_stalled_step is None else str(method.first_stalled_step)
            print(
                f"{method.method:<24}"
                f" {method.final_time:>12.6f}"
                f" {method.final_error:>14.3e}"
                f" {method.max_error:>13.3e}"
                f" {stalled:>20}"
                f" {method.unique_increment_count:>19}"
            )
        print()


def print_observer_guard_summary(results: tuple[ObserverGuardScenarioResult, ...]) -> None:
    print("Observer safeguard prototypes")
    print("scenario                              baseline                   events   guarded setting                 events")
    for result in results:
        print(
            f"{result.name:<35}"
            f" {result.baseline_label:<25}"
            f" {result.baseline_event_count:>6}"
            f"   {result.guard_label:<28}"
            f" {result.guarded_event_count:>6}"
        )
        print(f"  note: {result.note}")
    print()


def print_local_extremum_summary(result: LocalExtremumPrototypeResult) -> None:
    print("Local-extremum prototype")
    print(f"scenario: {result.name}")
    print(
        f" sample-argmax mean peak-time error = {result.argmax_mean_time_error:.3e}, max = {result.argmax_max_time_error:.3e}"
    )
    print(
        f" quadratic-vertex mean peak-time error = {result.quadratic_mean_time_error:.3e}, max = {result.quadratic_max_time_error:.3e}"
    )
    print(
        f" sample-argmax mean peak-value error = {result.argmax_mean_value_error:.3e}, quadratic-vertex = {result.quadratic_mean_value_error:.3e}"
    )
    print(f" note: {result.note}")
    print()


def print_threshold_interpolation_summary(result: ThresholdInterpolationPrototypeResult) -> None:
    print("Threshold-crossing interpolation prototype")
    print(f"scenario: {result.name}")
    print(
        f" sample-time mean crossing error = {result.sample_mean_time_error:.3e}, max = {result.sample_max_time_error:.3e}"
    )
    print(
        f" linear-inverse mean crossing error = {result.linear_mean_time_error:.3e}, max = {result.linear_max_time_error:.3e}"
    )
    print(
        f" Hermite mean crossing error = {result.hermite_mean_time_error:.3e}, max = {result.hermite_max_time_error:.3e}"
    )
    print(f" note: {result.note}")
    print()


def plot_results(
    mean_results: tuple[MeanScenarioResult, ...],
    time_results: tuple[TimeScenarioResult, ...],
    observer_results: tuple[ObserverGuardScenarioResult, ...],
    threshold_interpolation_result: ThresholdInterpolationPrototypeResult,
    local_extremum_result: LocalExtremumPrototypeResult,
) -> None:
    if plt is None:
        return

    figure, axes = plt.subplots(3, 2, figsize=(13, 13))
    axes = axes.ravel()

    labels = [result.name for result in mean_results]
    x = np.arange(len(labels))
    width = 0.35
    axes[0].bar(
        x - width / 2,
        [abs(result.running_error) for result in mean_results],
        width=width,
        label="runningMeanTime",
    )
    axes[0].bar(
        x + width / 2,
        [abs(result.compensated_error) for result in mean_results],
        width=width,
        label="compensated integral",
    )
    axes[0].set_xticks(x, labels, rotation=12, ha="right")
    axes[0].set_yscale("log")
    axes[0].set_ylabel("absolute error")
    axes[0].set_title("Time-weighted mean error in float32")
    axes[0].legend()

    zero_origin = next(result for result in time_results if result.name == "zero-origin drift")
    step_axis = np.arange(zero_origin.n_steps + 1, dtype=np.float64)
    reference = zero_origin.t0 + step_axis * zero_origin.dt
    traces = (
        ("direct add", direct_time_updates_float32(zero_origin.t0, zero_origin.dt, zero_origin.n_steps)),
        ("Kahan add", kahan_time_updates_float32(zero_origin.t0, zero_origin.dt, zero_origin.n_steps)),
        (
            "structured t0 + step*dt",
            structured_time_updates_float32(zero_origin.t0, zero_origin.dt, zero_origin.n_steps),
        ),
    )
    for label, trace in traces:
        axes[1].plot(step_axis, trace.astype(np.float64) - reference, label=label)
    axes[1].set_xlabel("step")
    axes[1].set_ylabel("time error")
    axes[1].set_title("Float32 time accumulation error")
    axes[1].legend()

    observer_labels = [result.name for result in observer_results]
    observer_x = np.arange(len(observer_labels))
    axes[2].bar(
        observer_x - width / 2,
        [result.baseline_event_count for result in observer_results],
        width=width,
        label="unguarded",
    )
    axes[2].bar(
        observer_x + width / 2,
        [result.guarded_event_count for result in observer_results],
        width=width,
        label="guarded",
    )
    axes[2].set_xticks(observer_x, observer_labels, rotation=12, ha="right")
    axes[2].set_ylabel("counted up-events")
    axes[2].set_title("Observer safeguards on synthetic traces")
    axes[2].legend()

    local_labels = ["time error", "value error"]
    local_x = np.arange(len(local_labels))
    axes[3].bar(
        np.arange(3) - width / 2,
        [
            threshold_interpolation_result.sample_mean_time_error,
            threshold_interpolation_result.linear_mean_time_error,
            threshold_interpolation_result.hermite_mean_time_error,
        ],
        width=width,
        color=["#c44e52", "#4c72b0", "#55a868"],
    )
    axes[3].set_xticks(np.arange(3), ["sample", "linear", "Hermite"])
    axes[3].set_yscale("log")
    axes[3].set_ylabel("mean absolute time error")
    axes[3].set_title("Threshold timestamp interpolation")

    axes[4].bar(
        local_x - width / 2,
        [local_extremum_result.argmax_mean_time_error, local_extremum_result.argmax_mean_value_error],
        width=width,
        label="sample argmax",
    )
    axes[4].bar(
        local_x + width / 2,
        [
            local_extremum_result.quadratic_mean_time_error,
            local_extremum_result.quadratic_mean_value_error,
        ],
        width=width,
        label="quadratic vertex",
    )
    axes[4].set_xticks(local_x, local_labels)
    axes[4].set_yscale("log")
    axes[4].set_ylabel("mean absolute error")
    axes[4].set_title("Three-sample local-max prototype")
    axes[4].legend()

    axes[5].axis("off")

    figure.tight_layout()
    if plt.get_backend().lower() != "agg":
        plt.show()


def main() -> None:
    mean_results = build_mean_scenarios()
    mean_origin_results = build_mean_origin_scenarios()
    time_results = build_time_scenarios()
    observer_results = build_observer_guard_scenarios()
    threshold_interpolation_result = build_threshold_interpolation_prototype()
    local_extremum_result = build_local_extremum_prototype()

    print_mean_summary(mean_results)
    print_mean_origin_summary(mean_origin_results)
    print_time_summary(time_results)
    print_observer_guard_summary(observer_results)
    print_threshold_interpolation_summary(threshold_interpolation_result)
    print_local_extremum_summary(local_extremum_result)
    print("Interpretation")
    print("- runningMeanTime can lose late small contributions once its incremental update drops below float32 resolution.")
    print("- compensated integral accumulation preserves those contributions much better because it sums at the scale of dt * value and divides once at the end.")
    print("- for time-weighted features, a large absolute t0 can break the elapsed-time denominator when it is formed as ti - t_start in float32; relative elapsed bookkeeping is much safer.")
    print("- direct time updates accumulate rounding error from the current absolute time; writing t = t + dt instead of t += dt does not change that.")
    print("- Kahan-style and structured time updates improve long-horizon endpoint accuracy, but once dt is below ulp(t) they still move in quantized jumps rather than resolving every step.")
    print("- reconstructing time from a step counter avoids repeated-add drift for fixed-step methods, but it is still limited by float32 spacing and, after 2^24, by float32 step-counter casts unless the counter stays wider than float32.")
    print("- for oscillation-oriented observers, amplitude floors, Schmitt-trigger-style hysteresis, and derivative thresholds are practical safeguards against overcounting tiny, chatter-driven, or noisy crossings.")
    print("- threshold timestamps have the same kind of interpolation tradeoff: sample times are coarse, inverse-linear interpolation is cheap and much better, and slope-aware Hermite interpolation is promising but needs stronger robustness checks before broad kernel use.")
    print("- for local-extremum observers, smaller dt or double precision are the safe user-side choices today when accurate peak timing matters; a shared interpolation helper is a strong candidate for future improvement.")

    plot_results(
        mean_results,
        time_results,
        observer_results,
        threshold_interpolation_result,
        local_extremum_result,
    )


if __name__ == "__main__":
    main()