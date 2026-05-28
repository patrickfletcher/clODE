import matplotlib.pyplot as plt
import numpy as np


def refine_quadratic_vertex(times: np.ndarray, values: np.ndarray) -> tuple[float, float]:
    coefficients = np.polyfit(times, values, deg=2)
    a, b, c = coefficients
    if np.isclose(a, 0.0):
        center_index = len(times) // 2
        return float(times[center_index]), float(values[center_index])
    vertex_time = float(np.clip(-b / (2.0 * a), times[0], times[-1]))
    vertex_value = float(a * vertex_time * vertex_time + b * vertex_time + c)
    return vertex_time, vertex_value


def local_extrema(
    sample_times: np.ndarray,
    sample_values: np.ndarray,
    sample_slopes: np.ndarray,
):
    maxima: list[tuple[float, float]] = []
    minima: list[tuple[float, float]] = []

    for index in range(2, len(sample_times)):
        window_times = sample_times[index - 2 : index + 1]
        window_values = sample_values[index - 2 : index + 1]
        slope_left = sample_slopes[index - 1]
        slope_right = sample_slopes[index]
        if slope_left > 0.0 and slope_right < 0.0:
            maxima.append(refine_quadratic_vertex(window_times, window_values))
        if slope_left < 0.0 and slope_right > 0.0:
            minima.append(refine_quadratic_vertex(window_times, window_values))

    return maxima, minima


def main() -> None:
    sample_times = np.linspace(0.0, 4.0 * np.pi, 29)
    sample_values = np.sin(sample_times)
    sample_slopes = np.cos(sample_times)
    dense_times = np.linspace(sample_times[0], sample_times[-1], 2000)
    dense_values = np.sin(dense_times)
    dense_slopes = np.cos(dense_times)

    maxima, minima = local_extrema(sample_times, sample_values, sample_slopes)

    figure, (axis_value, axis_slope) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    axis_value.plot(dense_times, dense_values, color="black", label="dense trace")
    axis_value.plot(sample_times, sample_values, "o", color="0.45", label="sampled trace")
    axis_value.plot(
        [time for time, _ in maxima],
        [value for _, value in maxima],
        "v",
        color="tab:red",
        label="refined maxima",
    )
    axis_value.plot(
        [time for time, _ in minima],
        [value for _, value in minima],
        "^",
        color="tab:blue",
        label="refined minima",
    )
    axis_value.set_ylabel("x")
    axis_value.set_title(
        "local_extremum triggers on sampled slope sign changes and refines with a three-sample quadratic fit"
    )
    axis_value.legend(loc="upper right")

    axis_slope.plot(dense_times, dense_slopes, color="black", label="dense dx/dt")
    axis_slope.plot(sample_times, sample_slopes, "o", color="0.45", label="sampled dx/dt")
    axis_slope.axhline(0.0, color="0.7", linestyle="--")
    axis_slope.set_xlabel("time")
    axis_slope.set_ylabel("dx/dt")
    axis_slope.legend(loc="upper right")

    figure.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
