from __future__ import annotations

from pathlib import Path

import numpy as np

import clode

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL = REPO_ROOT / "test/core_numerics/models/stable_linear.cl"

T_SPAN = (0.0, 1.0)
DT = 0.01
MAX_STEPS = max(32, int((T_SPAN[1] - T_SPAN[0]) / DT) + 1)

print(MAX_STEPS)

# Use the NVIDIA GPU:
PLATFORM = 1
DEVICE = 0

def main() -> None:
    simulator = clode.Simulator(
        src_file=str(MODEL),
        variables={"x": 2.0, "y": -1.5},
        parameters={"a": 0.5, "b": 1.25},
        stepper=clode.Stepper.rk4,
        dt=DT,
        dtmax=DT,
        t_span=T_SPAN,
        max_steps=MAX_STEPS,
        platform_id=PLATFORM,
        device_id=DEVICE,
    )
    simulator.transient()
    final_time = float(simulator.get_final_time()[0])

    trajectory_simulator = clode.TrajectorySimulator(
        src_file=str(MODEL),
        variables={"x": 2.0, "y": -1.5},
        parameters={"a": 0.5, "b": 1.25},
        stepper=clode.Stepper.rk4,
        dt=DT,
        dtmax=DT,
        t_span=T_SPAN,
        max_steps=MAX_STEPS,
        max_store=MAX_STEPS,
        nout=1,
        platform_id=PLATFORM,
        device_id=DEVICE,
    )
    trajectory = trajectory_simulator.trajectory()
    trajectory_final_time = float(trajectory_simulator.get_final_time()[0])

    print("requested t_span end:", T_SPAN[1])
    print("transient final_time:", final_time)
    print("trajectory final_time:", trajectory_final_time)
    print("trajectory stored times:", np.asarray(trajectory.t, dtype=np.float64))

    if final_time <= T_SPAN[1] + 1e-8:
        raise SystemExit("Fixed-step endpoint overshoot did not reproduce")

    print("issue reproduced: fixed-step integration steps past the requested end time")


if __name__ == "__main__":
    main()
