from __future__ import annotations

from typing import List

import clode
import matplotlib.pyplot as plt
import numpy as np


WINDOW_DURATION = 0.5
TOTAL_DURATION = 1.0
DT = 0.125
ATOL = 1e-12
MAX_STEPS = 64
MAX_STORE = 64


def stable_linear(
    t: float,
    variables: List[float],
    parameters: List[float],
    derivatives: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
     x: float = variables[0]
     y: float = variables[1]
     a: float = parameters[0]
     b: float = parameters[1]
     derivatives[0] = -a * x
     derivatives[1] = -b * y

def make_trajectory_simulator(t_span: tuple[float, float]) -> clode.TrajectorySimulator:
    return clode.TrajectorySimulator(
        rhs_equation=stable_linear,
        variables={"x": 2.0, "y": -1.5},
        parameters={"a": 0.5, "b": 1.25},
        stepper=clode.Stepper.rk4,
        t_span=t_span,
        dt=DT,
        dtmax=DT,
        max_steps=MAX_STEPS,
        max_store=MAX_STORE,
        single_precision=False,
    )


def make_feature_simulator(t_span: tuple[float, float]) -> clode.FeatureSimulator:
    return clode.FeatureSimulator(
        rhs_equation=stable_linear,
        variables={"x": 2.0, "y": -1.5},
        parameters={"a": 0.5, "b": 1.25},
        observer=clode.Observer.summary,
        stepper=clode.Stepper.rk4,
        t_span=t_span,
        dt=DT,
        dtmax=DT,
        max_steps=MAX_STEPS,
        single_precision=False,
    )


def concatenate_trajectory_segments(
    first: clode.TrajectoryOutput,
    second: clode.TrajectoryOutput,
) -> tuple[np.ndarray, np.ndarray]:
    time = np.concatenate([first.t, second.t[1:]])
    states = np.concatenate(
        [first.to_ndarray("x"), second.to_ndarray("x")[1:]],
        axis=0,
    )
    return time, states


def main() -> None:
    full_trajectory_simulator = make_trajectory_simulator((0.0, TOTAL_DURATION))
    split_trajectory_simulator = make_trajectory_simulator((0.0, WINDOW_DURATION))

    full_trajectory = full_trajectory_simulator.trajectory()
    first_window = split_trajectory_simulator.trajectory()
    split_trajectory_simulator.advance_tspan_to_attained_final_time()
    second_window = split_trajectory_simulator.trajectory()

    stitched_time, stitched_state = concatenate_trajectory_segments(
        first_window,
        second_window,
    )

    np.testing.assert_allclose(
        stitched_time,
        full_trajectory.t,
        atol=ATOL,
        rtol=0.0,
    )
    np.testing.assert_allclose(
        stitched_state,
        full_trajectory.to_ndarray("x"),
        atol=ATOL,
        rtol=0.0,
    )

    full_feature_simulator = make_feature_simulator((0.0, TOTAL_DURATION))
    split_feature_simulator = make_feature_simulator((0.0, WINDOW_DURATION))

    full_features = full_feature_simulator.features()
    split_feature_simulator.features()
    split_feature_simulator.advance_tspan_to_attained_final_time()
    split_features = split_feature_simulator.features()

    assert full_features is not None
    assert split_features is not None

    np.testing.assert_allclose(
        split_features.to_ndarray(),
        full_features.to_ndarray(),
        atol=ATOL,
        rtol=0.0,
    )

    trajectory_difference = np.max(
        np.abs(stitched_state - full_trajectory.to_ndarray("x"))
    )
    feature_difference = np.max(
        np.abs(split_features.to_ndarray() - full_features.to_ndarray())
    )

    print(f"max trajectory difference: {trajectory_difference:.3e}")
    print(f"max feature difference: {feature_difference:.3e}")
    print(f"mean x feature: {float(full_features.get_var_mean('x')):.12f}")
    print(f"mean y feature: {float(full_features.get_var_mean('y')):.12f}")

    figure, axis = plt.subplots(figsize=(8, 4.5))
    axis.plot(full_trajectory.t, full_trajectory.x["x"], label="single window", linewidth=2)
    axis.plot(
        stitched_time,
        stitched_state[:, 0],
        "--",
        label="split continuation",
        linewidth=2,
    )
    axis.set_xlabel("time")
    axis.set_ylabel("x")
    axis.set_title("Trajectory continuation matches a single long run")
    axis.legend()
    figure.tight_layout()
    if plt.get_backend().lower() != "agg":
        plt.show()


if __name__ == "__main__":
    main()