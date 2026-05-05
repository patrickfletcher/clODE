from __future__ import annotations

from pathlib import Path

import numpy as np

import clode

REPO_ROOT = Path(__file__).resolve().parents[1]
MODEL = REPO_ROOT / "test/core_numerics/models/stable_linear_aux.cl"


def feature_matrix(output: clode.ObserverOutput) -> np.ndarray:
    return np.asarray(output.to_ndarray(), dtype=np.float64)


def main() -> None:
    full = clode.FeatureSimulator(
        src_file=str(MODEL),
        variables={"x": 2.0, "y": -1.5},
        parameters={"a": 0.5, "b": 1.25},
        aux=["sum", "combo"],
        observer=clode.Observer.basic_all_variables,
        stepper=clode.Stepper.rk4,
        dt=0.125,
        dtmax=0.125,
        t_span=(0.0, 1.0),
        max_steps=64,
    )
    full_output = full.features()

    split = clode.FeatureSimulator(
        src_file=str(MODEL),
        variables={"x": 2.0, "y": -1.5},
        parameters={"a": 0.5, "b": 1.25},
        aux=["sum", "combo"],
        observer=clode.Observer.basic_all_variables,
        stepper=clode.Stepper.rk4,
        dt=0.125,
        dtmax=0.125,
        t_span=(0.0, 0.5),
        max_steps=64,
    )
    split.features()
    first_final_time = float(split.get_final_time()[0])
    split.set_tspan((first_final_time, 1.0))
    split_output = split.features()

    assert full_output is not None
    assert split_output is not None

    feature_names = full_output.get_feature_names()
    full_matrix = feature_matrix(full_output).reshape(-1)
    split_matrix = feature_matrix(split_output).reshape(-1)
    diff = np.abs(split_matrix - full_matrix)
    order = np.argsort(diff)[::-1]

    print("full final_time:", float(full.get_final_time()[0]))
    print("split first final_time:", first_final_time)
    print("split final_time after continuation:", float(split.get_final_time()[0]))
    print("max feature difference:", float(diff.max()))
    print("largest discrepancies:")
    for index in order[:5]:
        print(
            f"  {feature_names[index]}: split={split_matrix[index]:.9f}, "
            f"full={full_matrix[index]:.9f}, diff={diff[index]:.9f}"
        )

    if float(diff.max()) < 1e-3:
        raise SystemExit("basicall continuation mismatch did not reproduce")

    print("issue reproduced: basicall split-window continuation does not match a single long run")


if __name__ == "__main__":
    main()
