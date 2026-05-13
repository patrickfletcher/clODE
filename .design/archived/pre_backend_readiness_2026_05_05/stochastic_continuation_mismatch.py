from __future__ import annotations

from pathlib import Path

import numpy as np

import clode

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL = REPO_ROOT / "test/core_numerics/models/ornstein_uhlenbeck.cl"


def make_simulator(t_span: tuple[float, float]) -> clode.Simulator:
    simulator = clode.Simulator(
        src_file=str(MODEL),
        variables={"x": 0.0},
        parameters={"mu": 1.0, "sigma": 0.5},
        num_noise=1,
        stepper=clode.Stepper.stochastic_euler,
        dt=0.125,
        dtmax=0.125,
        t_span=t_span,
        max_steps=64,
    )
    simulator.set_repeat_ensemble(128)
    simulator.seed_rng(321)
    return simulator


def main() -> None:
    full = make_simulator((0.0, 1.0))
    split = make_simulator((0.0, 0.5))

    full.transient()

    split.transient()
    first_final_time = float(split.get_final_time().reshape(-1)[0])
    split.set_tspan((first_final_time, 1.0))
    split.transient()

    full_state = full.get_final_state().reshape(-1)
    split_state = split.get_final_state().reshape(-1)
    diff = np.abs(split_state - full_state)

    print("full final_time[0]:", float(full.get_final_time().reshape(-1)[0]))
    print("split first final_time[0]:", first_final_time)
    print("split final_time[0] after continuation:", float(split.get_final_time().reshape(-1)[0]))
    print("max state difference:", float(diff.max()))
    print("mean absolute difference:", float(diff.mean()))
    print("first five full states:", full_state[:5])
    print("first five split states:", split_state[:5])

    if float(diff.max()) < 1e-4:
        raise SystemExit("stochastic continuation mismatch did not reproduce")

    print("issue reproduced: split stochastic continuation does not match a single seeded run")


if __name__ == "__main__":
    main()
