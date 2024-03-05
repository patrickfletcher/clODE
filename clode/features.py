from __future__ import annotations

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from .function_converter import OpenCLRhsEquation
from .runtime import (
    CLDeviceType,
    CLVendor,
    FeatureSimulatorBase,
    ObserverParams,
    _clode_root_dir,
)
from .solver import Simulator, Stepper


class Observer(Enum):
    basic = "basic"
    basic_all_variables = "basicall"
    local_max = "localmax"
    neighbourhood_1 = "nhood1"
    neighbourhood_2 = "nhood2"
    threshold_2 = "thresh2"


# may want other functionality here.
# - full feature matrix, with names (pandas?)
# - separate diagnostics (period count, step count, min dt, max dt , ...), summarize on print along with observer info/pars?
class ObserverOutput:
    def __init__(
        self,
        observer_params: ObserverParams,
        results_array: np.ndarray[Any, np.dtype[np.float64]],
        num_result_features: int,
        variables: List[str],
        observer_type: Observer,
        feature_names: List[str],
    ):
        # print(type(observer_params))
        self._op = observer_params
        self._data = results_array
        self._num_result_features = num_result_features
        self._vars = variables
        self._observer_type = observer_type
        self._feature_names = feature_names

        shape = (
            self._num_result_features,
            len(results_array) // self._num_result_features,
        )
        self._data = results_array.reshape(shape).transpose()

    def get_feature_names(self) -> List[str]:
        return self._feature_names

    def _get_var(self, var: str, op: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        try:
            index = self._feature_names.index(f"{op} {var}")
            return self._data[:, index : index + 1].squeeze()
        except ValueError:
            raise NotImplementedError(
                f"{self._observer_type} does not track {op} {var}!"
            )

    def get_var_max(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(var, "max")

    def get_var_min(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(var, "min")

    def get_var_mean(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(var, "mean")

    def get_var_max_dt(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_max(f"d{var}/dt")

    def get_var_min_dt(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_min(f"d{var}/dt")

    def get_var_count(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var("count", var)


class FeatureSimulator(Simulator):
    """
    A class for computing features from a CLODE model.

    Parameters
    ----------
    src_file : str
        Path to the CLODE model source file.
    variable_names : List[str]
        List of variable names in the model.
    parameter_names : List[str]
        List of parameter names in the model.
    aux : List[str], optional
        List of auxiliary variable names in the model, by default None
    num_noise : int, optional
        Number of noise variables in the model, by default 1
    event_var : str, optional
        Name of the variable to use for event detection, by default ""
    feature_var : str, optional
        Name of the variable to use for feature detection, by default ""
    observer_max_event_count : int, optional
        Maximum number of events to detect, by default 100
    observer_min_x_amp : float, optional
        Minimum amplitude of the feature variable to detect, by default 1.0
    observer_min_imi : float, optional
        Minimum inter-event interval to detect, by default 1
    observer_neighbourhood_radius : float, optional
        Radius of the neighbourhood to use for event detection, by default 0.01
    observer_x_up_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        upper threshold, by default 0.3
    observer_x_down_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 0.2
    observer_dx_up_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        upper threshold, by default 0
    observer_dx_down_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 0
    observer_eps_dx : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 1e-7
    tspan : tuple[float, float], optional
        Time span for the simulation, by default (0.0, 1000.0)
    stepper : Stepper, optional
        Stepper to use for the simulation, by default Stepper.euler
    single_precision : bool, optional
        Whether to use single precision for the simulation, by default False
    dt : float, optional
        Time step for the simulation, by default 0.1
    dtmax : float, optional
        Maximum time step for the simulation, by default 1.0
    atol : float, optional
        Absolute tolerance for the simulation, by default 1e-6
    rtol : float, optional
        Relative tolerance for the simulation, by default 1e-6
    max_steps : int, optional
        Maximum number of steps for the simulation, by default 100000
    max_error : float, optional
        Maximum error for the simulation, by default 1e-3
    max_num_events : int, optional
        Maximum number of events to detect, by default 100
    min_x_amp : float, optional
        Minimum amplitude of the feature variable to detect, by default 1.0
    min_imi : float, optional
        Minimum inter-event interval to detect, by default 1
    neighbourhood_radius : float, optional
        Radius of the neighbourhood to use for event detection, by default 0.01
    x_up_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        upper threshold, by default 0.3
    x_down_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 0.2
    dx_up_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        upper threshold, by default 0
    dx_down_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 0
    eps_dx : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 1e-7
    max_event_count : int, optional
        Maximum number of events to detect, by default 100
    min_x_amp : float, optional
        Minimum amplitude of the feature variable to detect, by default 1.0
    min_imi : float, optional
        Minimum inter-event interval to detect, by default 1
    neighbourhood_radius : float, optional
        Radius of the neighbourhood to use for event detection, by default 0.01
    x_up_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        upper threshold, by default 0.3
    x_down_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 0.2
    dx_up_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        upper threshold, by default 0
    dx_down_thresh : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 0
    eps_dx : float, optional
        Threshold for detecting an event when the feature variable crosses the
        lower threshold, by default 1e-7

    Returns:
    --------
    CLODEFeatures
        A CLODEFeatures object.

    Examples
    --------
    >>> import clode
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> model = clode.FeatureSimulator(
    ...     src_file="examples/lorenz96.c",
    ...     variable_names=["x"],
    ...     parameter_names=["F"],

    ... )
    >>> model.set_parameter_values({"F": 8.0})
    >>> model.set_initial_values({"x": np.random.rand(40)})
    >>> model.simulate()
    >>> model.plot()
    >>> plt.show()"""

    _integrator: FeatureSimulatorBase

    def __init__(
        self,
        variables: Dict[str, float],
        parameters: Dict[str, float],
        src_file: Optional[str] = None,
        rhs_equation: Optional[OpenCLRhsEquation] = None,
        supplementary_equations: Optional[List[Callable[[Any], Any]]] = None,
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        t_span: Tuple[float, float] = (0.0, 1000.0),
        stepper: Stepper = Stepper.rk4,
        observer: Observer = Observer.basic_all_variables,
        single_precision: bool = True,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 10000000,
        max_store: int = 10000000,
        nout: int = 1,
        device_type: CLDeviceType | None = None,
        vendor: CLVendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
        device_ids: List[int] | None = None,
        event_var: str = "",
        feature_var: str = "",
        observer_max_event_count: int = 100,
        observer_min_x_amp: float = 0.0,
        observer_min_imi: float = 0.0,
        observer_neighbourhood_radius: float = 0.05,
        observer_x_up_thresh: float = 0.3,
        observer_x_down_thresh: float = 0.2,
        observer_dx_up_thresh: float = 0,
        observer_dx_down_thresh: float = 0,
        observer_eps_dx: float = 0.0,
    ) -> None:

        self._observer_type = observer

        event_var_idx = (
            list(variables.keys()).index(event_var) if event_var != "" else 0
        )
        feature_var_idx = (
            list(variables.keys()).index(feature_var) if feature_var != "" else 0
        )

        self._op = ObserverParams(
            event_var_idx,
            feature_var_idx,
            observer_max_event_count,
            observer_min_x_amp,
            observer_min_imi,
            observer_neighbourhood_radius,
            observer_x_up_thresh,
            observer_x_down_thresh,
            observer_dx_up_thresh,
            observer_dx_down_thresh,
            observer_eps_dx,
        )

        super().__init__(
            variables=variables,
            parameters=parameters,
            src_file=src_file,
            rhs_equation=rhs_equation,
            supplementary_equations=supplementary_equations,
            aux=aux,
            num_noise=num_noise,
            t_span=t_span,
            stepper=stepper,
            single_precision=single_precision,
            dt=dt,
            dtmax=dtmax,
            abstol=abstol,
            reltol=reltol,
            max_steps=max_steps,
            max_store=max_store,
            nout=nout,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )

    def _build_integrator(self) -> None:
        self._integrator = FeatureSimulatorBase(
            self._pi,
            self._stepper.value,
            self._observer_type.value,
            self._single_precision,
            self._runtime,
            _clode_root_dir,
        )

    def _init_integrator(self) -> None:
        vars_array, pars_array = self._pack_data()
        self._integrator.initialize(
            self.tspan,
            vars_array,
            pars_array,
            self._sp,
            self._op,
        )

    def get_feature_names(self) -> List[str]:
        return self._integrator.get_feature_names()

    def features(
        self, initialize_observer: Optional[bool] = None, update_x0: bool = True
    ) -> ObserverOutput:
        """Run a simulation with feature detection.

        Returns:
            None
        """
        if not self.is_initialized:
            raise RuntimeError("Simulator is not initialized")

        if initialize_observer is not None:
            print("Reinitializing observer")
            self._integrator.features(initialize_observer)
        else:
            self._integrator.features()
        if update_x0:
            self.shift_x0()

        return self.get_observer_results()

    def get_observer_results(self) -> ObserverOutput:
        """Get the feature data.

        Returns:
            np.array: The feature data.
        """
        self._result_integrator = self._integrator.get_f()
        self._num_result_features = self._integrator.get_n_features()
        return ObserverOutput(
            self._op,
            np.array(self._result_integrator),
            self._num_result_features,
            self.variable_names,
            self._observer_type,
            self._integrator.get_feature_names(),
        )
