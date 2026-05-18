from __future__ import annotations

import numpy as np

import clode
from test.core_numerics.helpers import (
    STABLE_LINEAR_PARAMETERS,
    STABLE_LINEAR_VARIABLES,
    device_kwargs_for_tests,
    make_simulator,
    model_path,
)


def test_simulator_accepts_ivp_and_matches_legacy_construction() -> None:
    ivp = clode.InitialValueProblem(
        src_file=model_path("stable_linear.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
    )
    legacy = make_simulator(
        "stable_linear",
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.05,
        max_steps=512,
    )
    via_ivp = clode.Simulator(
        ivp=ivp,
        stepper=clode.Stepper.rk4,
        t_span=(0.0, 1.0),
        dt=0.05,
        max_steps=512,
        single_precision=True,
        **device_kwargs_for_tests(),
    )

    legacy.transient(update_x0=False, fetch_results=False)
    via_ivp.transient(update_x0=False, fetch_results=False)

    assert via_ivp.ivp is ivp
    np.testing.assert_allclose(via_ivp.get_final_state(), legacy.get_final_state())