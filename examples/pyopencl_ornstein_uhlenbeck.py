from __future__ import annotations

import argparse
import os
from typing import List

os.environ["_CLODE_BACKEND"] = "pyopencl"

import clode
import numpy as np


def ornstein_uhlenbeck(
    t: float,
    variables: List[float],
    parameters: List[float],
    derivatives: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    x: float = variables[0]
    mu: float = parameters[0]
    sigma: float = parameters[1]
    dw: float = wiener[0]
    derivatives[0] = mu - x + sigma * dw


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Run an Ornstein-Uhlenbeck ensemble through the internal PyOpenCL "
            "backend selector used during the backend migration."
        )
    )
    parser.add_argument("--platform-id", type=int, default=None)
    parser.add_argument("--device-id", type=int, default=None)
    parser.add_argument("--ensemble-size", type=int, default=8192)
    parser.add_argument("--t-final", type=float, default=20.0)
    parser.add_argument("--dt", type=float, default=0.01)
    parser.add_argument("--mu", type=float, default=1.0)
    parser.add_argument("--sigma", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=1234)
    parser.add_argument(
        "--print-opencl",
        action="store_true",
        help="Print the OpenCL platforms visible to clODE before running.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if (args.platform_id is None) != (args.device_id is None):
        raise ValueError("Specify both --platform-id and --device-id together")

    if args.print_opencl:
        clode.print_opencl()

    integrator = clode.Simulator(
        rhs_equation=ornstein_uhlenbeck,
        variables={"x": 0.0},
        parameters={"mu": args.mu, "sigma": args.sigma},
        num_noise=1,
        stepper=clode.Stepper.stochastic_euler,
        single_precision=True,
        t_span=(0.0, args.t_final),
        dt=args.dt,
        platform_id=args.platform_id,
        device_id=args.device_id,
    )

    integrator.set_repeat_ensemble(args.ensemble_size)
    integrator.seed_rng(args.seed)
    integrator.transient()

    final_state = np.asarray(integrator.get_final_state(), dtype=np.float64).reshape(-1)
    expected_mean = args.mu
    expected_variance = args.sigma**2 / 2.0

    runtime_text = "auto-selected runtime"
    if args.platform_id is not None and args.device_id is not None:
        runtime_text = (
            f"platform_id={args.platform_id}, device_id={args.device_id}"
        )

    print(f"Backend selector: {os.environ['_CLODE_BACKEND']}")
    print(f"Runtime selection: {runtime_text}")
    print(f"Ensemble size: {args.ensemble_size}")
    print(f"Expected stationary mean: {expected_mean:.6f}")
    print(f"Simulated mean: {np.mean(final_state):.6f}")
    print(f"Expected stationary variance: {expected_variance:.6f}")
    print(f"Simulated variance: {np.var(final_state):.6f}")


if __name__ == "__main__":
    main()