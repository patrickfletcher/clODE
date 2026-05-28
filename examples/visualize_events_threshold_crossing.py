import matplotlib.pyplot as plt
import numpy as np


def crossing_times(times: np.ndarray, values: np.ndarray, threshold: float):
    rising: list[float] = []
    falling: list[float] = []
    for index in range(1, len(times)):
        t0 = times[index - 1]
        t1 = times[index]
        x0 = values[index - 1]
        x1 = values[index]
        x_interp = float(t0 + (threshold - x0) * (t1 - t0) / (x1 - x0))
        if x0 <= threshold and x1 > threshold:
            rising.append(x_interp)
        if x0 >= threshold and x1 < threshold:
            falling.append(x_interp)
    return rising, falling


def main() -> None:
    sample_times = np.linspace(0.0, 4.0 * np.pi, 33)
    sample_values = np.sin(sample_times)
    dense_times = np.linspace(sample_times[0], sample_times[-1], 2000)
    dense_values = np.sin(dense_times)

    normalized_threshold = 0.675
    threshold = sample_values.min() + normalized_threshold * (
        sample_values.max() - sample_values.min()
    )
    rising, falling = crossing_times(sample_times, sample_values, threshold)

    figure, axis = plt.subplots(figsize=(10, 5))
    axis.plot(dense_times, dense_values, color="black", label="dense trace")
    axis.plot(sample_times, sample_values, "o", color="0.45", label="sampled trace")
    axis.axhline(threshold, color="tab:green", linestyle="--", label="threshold")

    for index, time in enumerate(rising):
        label = "rising event" if index == 0 else None
        axis.axvline(time, color="tab:red", linestyle=":", label=label)
    for index, time in enumerate(falling):
        label = "falling event" if index == 0 else None
        axis.axvline(time, color="tab:blue", linestyle=":", label=label)

    axis.set_title(
        "threshold_crossing and normalized_threshold_crossing share the same crossing rule"
    )
    axis.set_xlabel("time")
    axis.set_ylabel("x")
    axis.text(
        0.02,
        0.03,
        f"normalized threshold = {normalized_threshold:.3f}\nabsolute threshold = {threshold:.3f}",
        transform=axis.transAxes,
        ha="left",
        va="bottom",
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )
    axis.legend(loc="upper right")
    figure.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
