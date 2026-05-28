import matplotlib.pyplot as plt
import numpy as np


def schmitt_transitions(
    times: np.ndarray,
    values: np.ndarray,
    x_up: float,
    x_down: float,
):
    up_times: list[float] = []
    down_times: list[float] = []
    in_upstate = bool(values[0] > x_up)

    for index in range(1, len(times)):
        t0 = times[index - 1]
        t1 = times[index]
        x0 = values[index - 1]
        x1 = values[index]
        if not in_upstate and x0 <= x_up and x1 > x_up:
            up_times.append(float(np.interp(x_up, [x0, x1], [t0, t1])))
            in_upstate = True
        elif in_upstate and x0 >= x_down and x1 < x_down:
            down_times.append(float(np.interp(x_down, [x0, x1], [t0, t1])))
            in_upstate = False

    return up_times, down_times


def main() -> None:
    sample_times = np.linspace(0.0, 6.0 * np.pi, 61)
    sample_values = np.sin(sample_times)
    dense_times = np.linspace(sample_times[0], sample_times[-1], 3000)
    dense_values = np.sin(dense_times)

    x_up = 0.45
    x_down = -0.15
    up_times, down_times = schmitt_transitions(sample_times, sample_values, x_up, x_down)

    figure, axis = plt.subplots(figsize=(10, 5))
    axis.plot(dense_times, dense_values, color="black", label="dense trace")
    axis.plot(sample_times, sample_values, "o", color="0.45", label="sampled trace")
    axis.axhline(x_up, color="tab:red", linestyle="--", label="x_up")
    axis.axhline(x_down, color="tab:blue", linestyle="--", label="x_down")

    for index, time in enumerate(up_times):
        label = "up transition" if index == 0 else None
        axis.axvline(time, color="tab:red", linestyle=":", label=label)
    for index, time in enumerate(down_times):
        label = "down transition" if index == 0 else None
        axis.axvline(time, color="tab:blue", linestyle=":", label=label)

    for up_time, down_time in zip(up_times, down_times):
        axis.axvspan(up_time, down_time, color="tab:red", alpha=0.12)

    axis.set_title(
        "schmitt_trigger uses an x-only two-boundary state machine; x_up == x_down is allowed"
    )
    axis.set_xlabel("time")
    axis.set_ylabel("x")
    axis.legend(loc="upper right")
    figure.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
