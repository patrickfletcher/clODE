from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import clode
import numpy as np

MODEL_DIR = Path(__file__).parent / "models"

STABLE_LINEAR_VARIABLES = {"x": 2.0, "y": -1.5}
STABLE_LINEAR_PARAMETERS = {"a": 0.5, "b": 1.25}
STABLE_LINEAR_AUX = ["sum", "combo"]

OU_VARIABLES = {"x": 0.0}
OU_PARAMETERS = {"mu": 1.0, "sigma": 0.5}

HOPF_VARIABLES = {"x": 1.0, "y": 0.0}
HOPF_PARAMETERS = {"mu": 1.0, "omega": 1.0}
HOPF_AUX = ["r2"]


def _optional_int_env(name: str) -> int | None:
    value = os.getenv(name)
    return None if value in (None, "") else int(value)


TEST_PLATFORM_ID = _optional_int_env("CLODE_TEST_PLATFORM_ID")
TEST_DEVICE_ID = _optional_int_env("CLODE_TEST_DEVICE_ID")


def model_path(name: str) -> str:
    return str(MODEL_DIR / name)


def _with_test_device(kwargs: dict[str, Any]) -> dict[str, Any]:
    options = dict(kwargs)
    if TEST_PLATFORM_ID is not None and "platform_id" not in options:
        options["platform_id"] = TEST_PLATFORM_ID
    if TEST_DEVICE_ID is not None and "device_id" not in options:
        options["device_id"] = TEST_DEVICE_ID
    return options


def device_kwargs_for_tests(**kwargs: Any) -> dict[str, Any]:
    return _with_test_device(kwargs)


def _problem(name: str) -> dict[str, Any]:
    if name == "stable_linear":
        return {
            "src_file": model_path("stable_linear.cl"),
            "variables": STABLE_LINEAR_VARIABLES.copy(),
            "parameters": STABLE_LINEAR_PARAMETERS.copy(),
            "aux": [],
            "num_noise": 0,
        }
    if name == "stable_linear_aux":
        return {
            "src_file": model_path("stable_linear_aux.cl"),
            "variables": STABLE_LINEAR_VARIABLES.copy(),
            "parameters": STABLE_LINEAR_PARAMETERS.copy(),
            "aux": STABLE_LINEAR_AUX.copy(),
            "num_noise": 0,
        }
    if name == "ornstein_uhlenbeck":
        return {
            "src_file": model_path("ornstein_uhlenbeck.cl"),
            "variables": OU_VARIABLES.copy(),
            "parameters": OU_PARAMETERS.copy(),
            "aux": [],
            "num_noise": 1,
        }
    if name == "hopf_normal_form":
        return {
            "src_file": model_path("hopf_normal_form.cl"),
            "variables": HOPF_VARIABLES.copy(),
            "parameters": HOPF_PARAMETERS.copy(),
            "aux": HOPF_AUX.copy(),
            "num_noise": 0,
        }
    raise ValueError(f"Unknown model: {name}")


def make_simulator(
    model_name: str,
    *,
    stepper: clode.Stepper,
    t_span: tuple[float, float],
    dt: float,
    dtmax: float | None = None,
    max_steps: int = 100000,
    single_precision: bool = True,
    **kwargs: Any,
) -> clode.Simulator:
    problem = _problem(model_name)
    options = _with_test_device(kwargs)
    return clode.Simulator(
        src_file=problem["src_file"],
        variables=problem["variables"],
        parameters=problem["parameters"],
        aux=problem["aux"],
        num_noise=problem["num_noise"],
        stepper=stepper,
        dt=dt,
        dtmax=dt if dtmax is None else dtmax,
        t_span=t_span,
        max_steps=max_steps,
        single_precision=single_precision,
        **options,
    )


def make_trajectory_simulator(
    model_name: str,
    *,
    stepper: clode.Stepper,
    t_span: tuple[float, float],
    dt: float,
    dtmax: float | None = None,
    max_steps: int = 100000,
    max_store: int = 1000,
    nout: int = 1,
    single_precision: bool = True,
    **kwargs: Any,
) -> clode.TrajectorySimulator:
    problem = _problem(model_name)
    options = _with_test_device(kwargs)
    return clode.TrajectorySimulator(
        src_file=problem["src_file"],
        variables=problem["variables"],
        parameters=problem["parameters"],
        aux=problem["aux"],
        num_noise=problem["num_noise"],
        stepper=stepper,
        dt=dt,
        dtmax=dt if dtmax is None else dtmax,
        t_span=t_span,
        max_steps=max_steps,
        max_store=max_store,
        nout=nout,
        single_precision=single_precision,
        **options,
    )


def make_feature_simulator(
    model_name: str,
    *,
    stepper: clode.Stepper,
    t_span: tuple[float, float],
    dt: float,
    dtmax: float | None = None,
    max_steps: int = 100000,
    single_precision: bool = True,
    **kwargs: Any,
) -> clode.FeatureSimulator:
    problem = _problem(model_name)
    options = _with_test_device(kwargs)
    return clode.FeatureSimulator(
        src_file=problem["src_file"],
        variables=problem["variables"],
        parameters=problem["parameters"],
        aux=problem["aux"],
        num_noise=problem["num_noise"],
        observer=clode.Observer.summary,
        stepper=stepper,
        dt=dt,
        dtmax=dt if dtmax is None else dtmax,
        t_span=t_span,
        max_steps=max_steps,
        single_precision=single_precision,
        **options,
    )


def set_repeat_ensemble_and_seed(
    simulator: clode.Simulator,
    *,
    ensemble_size: int,
    seed: int,
) -> None:
    simulator.set_repeat_ensemble(ensemble_size)
    simulator.seed_rng(seed)


def structured_to_array(output: Any, field: str) -> np.ndarray:
    if field == "x":
        return np.asarray(output.to_ndarray("x"), dtype=np.float64)
    if field == "dx":
        return np.asarray(output.to_ndarray("dx"), dtype=np.float64)
    if field == "aux":
        return np.asarray(output.to_ndarray("aux"), dtype=np.float64)
    raise ValueError(f"Unknown trajectory field: {field}")
