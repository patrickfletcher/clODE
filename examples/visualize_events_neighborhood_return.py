import math

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


def interpolate_ball_exit(
    t0: float,
    t1: float,
    start_point: tuple[float, float],
    end_point: tuple[float, float],
    anchor_point: tuple[float, float],
    x_range: float,
    y_range: float,
    radius: float,
) -> tuple[float, tuple[float, float]]:
    a = 0.0
    b = 0.0
    c = -(radius * radius)
    for start_value, end_value, anchor_value, axis_range in (
        (start_point[0], end_point[0], anchor_point[0], x_range),
        (start_point[1], end_point[1], anchor_point[1], y_range),
    ):
        if np.isclose(axis_range, 0.0):
            continue
        normalized_start = (start_value - anchor_value) / axis_range
        normalized_delta = (end_value - start_value) / axis_range
        a += normalized_delta * normalized_delta
        b += 2.0 * normalized_start * normalized_delta
        c += normalized_start * normalized_start

    alpha = 1.0
    if not np.isclose(t1, t0):
        if not np.isclose(a, 0.0):
            discriminant = max(b * b - 4.0 * a * c, 0.0)
            alpha = (-b + math.sqrt(discriminant)) / (2.0 * a)
        elif not np.isclose(b, 0.0):
            alpha = -c / b
    alpha = float(np.clip(alpha, 0.0, 1.0))

    event_time = t0 + alpha * (t1 - t0)
    event_point = (
        start_point[0] + alpha * (end_point[0] - start_point[0]),
        start_point[1] + alpha * (end_point[1] - start_point[1]),
    )
    return event_time, event_point


def neighborhood_return_events(
    times: np.ndarray,
    x_values: np.ndarray,
    y_values: np.ndarray,
    anchor_threshold: float,
    radius: float,
):
    x_min = float(np.min(x_values))
    x_max = float(np.max(x_values))
    y_min = float(np.min(y_values))
    y_max = float(np.max(y_values))
    x_range = x_max - x_min
    y_range = y_max - y_min
    x_threshold = x_min + anchor_threshold * x_range

    found_anchor = False
    anchor_index = -1
    anchor_point = (0.0, 0.0)
    in_neighborhood = False
    event_times: list[float] = []
    event_points: list[tuple[float, float]] = []

    for index in range(1, len(times)):
        if not found_anchor:
            if x_values[index - 1] > x_threshold and x_values[index] < x_threshold:
                found_anchor = True
                anchor_index = index
                anchor_point = (float(x_values[index]), float(y_values[index]))
                in_neighborhood = True
            continue

        normalized_dx = 0.0 if np.isclose(x_range, 0.0) else (x_values[index] - anchor_point[0]) / x_range
        normalized_dy = 0.0 if np.isclose(y_range, 0.0) else (y_values[index] - anchor_point[1]) / y_range
        last_in_neighborhood = in_neighborhood
        in_neighborhood = normalized_dx * normalized_dx + normalized_dy * normalized_dy < radius * radius
        if last_in_neighborhood and not in_neighborhood:
            event_time, event_point = interpolate_ball_exit(
                float(times[index - 1]),
                float(times[index]),
                (float(x_values[index - 1]), float(y_values[index - 1])),
                (float(x_values[index]), float(y_values[index])),
                anchor_point,
                x_range,
                y_range,
                radius,
            )
            event_times.append(event_time)
            event_points.append(event_point)

    return {
        "anchor_index": anchor_index,
        "anchor_point": anchor_point,
        "event_times": event_times,
        "event_points": event_points,
        "x_threshold": x_threshold,
        "x_range": x_range,
        "y_range": y_range,
    }


def main() -> None:
    sample_times = np.linspace(0.0, 6.0 * np.pi, 121)
    sample_x = 0.8 * np.cos(sample_times) + 0.15 * np.cos(2.0 * sample_times)
    sample_y = 0.6 * np.sin(sample_times)

    anchor_threshold = 0.18
    radius = 0.2
    result = neighborhood_return_events(
        sample_times,
        sample_x,
        sample_y,
        anchor_threshold,
        radius,
    )

    anchor_index = result["anchor_index"]
    anchor_x, anchor_y = result["anchor_point"]
    event_times = result["event_times"]
    event_points = result["event_points"]

    figure, (axis_time, axis_phase) = plt.subplots(2, 1, figsize=(10, 8))
    axis_time.plot(sample_times, sample_x, color="black", label="event variable")
    axis_time.axhline(
        result["x_threshold"],
        color="tab:green",
        linestyle="--",
        label="anchor threshold",
    )
    if anchor_index >= 0:
        axis_time.axvline(
            sample_times[anchor_index],
            color="tab:purple",
            linestyle=":",
            label="anchor latch",
        )
    for index, event_time in enumerate(event_times):
        label = "neighborhood exit" if index == 0 else None
        axis_time.axvline(event_time, color="tab:red", linestyle=":", label=label)
    axis_time.set_title(
        "neighborhood_return keeps a sampled anchor x0 and linearly refines the normalized-ball exit time"
    )
    axis_time.set_xlabel("time")
    axis_time.set_ylabel("event variable")
    axis_time.legend(loc="upper right")

    axis_phase.plot(sample_x, sample_y, color="black", label="trajectory")
    if anchor_index >= 0:
        axis_phase.plot(anchor_x, anchor_y, "o", color="tab:purple", label="anchor x0")
        axis_phase.add_patch(
            Ellipse(
                (anchor_x, anchor_y),
                width=2.0 * radius * result["x_range"],
                height=2.0 * radius * result["y_range"],
                facecolor="tab:red",
                alpha=0.12,
                edgecolor="tab:red",
            )
        )
    if event_points:
        axis_phase.plot(
            [point[0] for point in event_points],
            [point[1] for point in event_points],
            "x",
            color="tab:red",
            markersize=8,
            label="interpolated exits",
        )
    axis_phase.set_xlabel("x")
    axis_phase.set_ylabel("y")
    axis_phase.set_aspect("equal", adjustable="box")
    axis_phase.legend(loc="upper right")

    figure.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
