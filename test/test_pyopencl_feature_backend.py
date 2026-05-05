from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pyopencl")

import clode
from test.core_numerics.helpers import (
    HOPF_AUX,
    HOPF_PARAMETERS,
    HOPF_VARIABLES,
    STABLE_LINEAR_AUX,
    STABLE_LINEAR_PARAMETERS,
    STABLE_LINEAR_VARIABLES,
    device_kwargs_for_tests,
    model_path,
)


FIXED_DT = 0.05
FIXED_MAX_STEPS = 512


def _make_pyopencl_feature(monkeypatch: pytest.MonkeyPatch, **kwargs: object) -> clode.FeatureSimulator:
    monkeypatch.setenv("_CLODE_BACKEND", "pyopencl")
    return clode.FeatureSimulator(**device_kwargs_for_tests(**kwargs))


def _make_cpp_feature(**kwargs: object) -> clode.FeatureSimulator:
    return clode.FeatureSimulator(**device_kwargs_for_tests(**kwargs))


def test_pyopencl_feature_backend_matches_cpp_for_basicall_statistics(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("_CLODE_BACKEND", raising=False)
    cpp = _make_cpp_feature(
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.basic_all_variables,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )
    pyopencl_backend = _make_pyopencl_feature(
        monkeypatch,
        src_file=model_path("stable_linear_aux.cl"),
        variables=STABLE_LINEAR_VARIABLES.copy(),
        parameters=STABLE_LINEAR_PARAMETERS.copy(),
        aux=STABLE_LINEAR_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.basic_all_variables,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 1.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    cpp_result = cpp.features(update_x0=False)
    pyopencl_result = pyopencl_backend.features(update_x0=False)

    assert cpp_result is not None
    assert pyopencl_result is not None
    assert pyopencl_result.get_feature_names() == cpp_result.get_feature_names()
    np.testing.assert_allclose(
        pyopencl_result.to_ndarray(),
        cpp_result.to_ndarray(),
        atol=1e-6,
        rtol=0.0,
    )


def test_pyopencl_feature_backend_localmax_rebuild_updates_event_storage(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    pyopencl_backend = _make_pyopencl_feature(
        monkeypatch,
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.local_max,
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 8.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    pyopencl_backend.set_observer_parameters(max_event_timestamps=3)

    pyopencl_result = pyopencl_backend.features(update_x0=False)

    assert pyopencl_result is not None
    assert len(
        [
            name
            for name in pyopencl_result.get_feature_names()
            if name.startswith("localmax event time")
        ]
    ) == 3
    assert "-DN_STORE_EVENTS=3" in pyopencl_backend._integrator.get_program_string()
    assert pyopencl_result.to_ndarray().shape[-1] == pyopencl_backend._integrator.get_n_features()


def test_pyopencl_feature_backend_matches_cpp_for_two_pass_threshold_observer(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("_CLODE_BACKEND", raising=False)
    cpp = _make_cpp_feature(
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.threshold_2,
        event_var="x",
        feature_var="x",
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 8.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )
    pyopencl_backend = _make_pyopencl_feature(
        monkeypatch,
        src_file=model_path("hopf_normal_form.cl"),
        variables=HOPF_VARIABLES.copy(),
        parameters=HOPF_PARAMETERS.copy(),
        aux=HOPF_AUX.copy(),
        num_noise=0,
        observer=clode.Observer.threshold_2,
        event_var="x",
        feature_var="x",
        stepper=clode.Stepper.rk4,
        dt=FIXED_DT,
        dtmax=FIXED_DT,
        t_span=(0.0, 8.0),
        max_steps=FIXED_MAX_STEPS,
        single_precision=True,
    )

    cpp_result = cpp.features(update_x0=False)
    pyopencl_result = pyopencl_backend.features(update_x0=False)

    assert cpp_result is not None
    assert pyopencl_result is not None
    assert pyopencl_result.get_feature_names() == cpp_result.get_feature_names()
    np.testing.assert_allclose(
        pyopencl_result.to_ndarray(),
        cpp_result.to_ndarray(),
        atol=1e-5,
        rtol=0.0,
    )