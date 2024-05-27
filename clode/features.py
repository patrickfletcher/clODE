from __future__ import annotations

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
from numpy.lib import recfunctions as rfn

from clode.cpp.clode_cpp_wrapper import (
    FeatureSimulatorBase,
    ObserverParams,
    SolverParams,
)

from .function_converter import OpenCLRhsEquation
from .runtime import CLDeviceType, CLVendor, _clode_root_dir
from .solver import Simulator, Stepper

# TODO[API]: defaults for ObserverParams are here and in wrapper. Should be only ONE place globally.
# - Prefer the struct defaults? Create a default config?
# TODO[API]: Observer param subsets are specific to different observers. Should model each as a class


class Observer(Enum):
    basic = "basic"
    basic_all_variables = "basicall"
    local_max = "localmax"
    neighbourhood_1 = "nhood1"
    neighbourhood_2 = "nhood2"
    threshold_2 = "thresh2"


# may want other functionality here.
# - full feature matrix, with names [numpy sturctured arrays?]
# - separate diagnostic features (period count, step count, min dt, max dt , ...)
class ObserverOutput:
    def __init__(
        self,
        observer_params: ObserverParams,
        feature_array: np.ndarray[Any, np.dtype[np.float64]],
        num_features: int,
        variables: List[str],
        observer_type: Observer,
        feature_names: List[str],
        ensemble_shape: Tuple,
    ) -> None:
        self._op = observer_params
        self._num_features = num_features
        self._vars = variables
        self._observer_type = observer_type
        self._feature_names = feature_names
        self._ensemble_shape = ensemble_shape

        # support indexing F via feature name directly
        F_dtype = np.dtype(
            {"names": feature_names, "formats": [np.float64] * len(feature_names)}
        )
        self.F = rfn.unstructured_to_structured(feature_array, dtype=F_dtype)

    def __repr__(self) -> str:
        ensemble_size = len(self.F[self._feature_names[0]])
        num_features = len(self._feature_names)
        feature_names = self._feature_names
        return f"ObserverOutput( ensemble size: {ensemble_size}, number of features: {num_features}, feature_names: {feature_names})"

    def to_ndarray(self, **kwargs):
        return rfn.structured_to_unstructured(self.F, **kwargs)

    def get_feature_names(self) -> List[str]:
        return self._feature_names

    def _get_var(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        try:
            result = self.F[var].squeeze().reshape(self._ensemble_shape)
            result = result[0] if result.size == 1 else result
            return result
        except ValueError:
            raise NotImplementedError(f"{self._observer_type} does not track {var}!")

    def get_var_max(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(" ".join(["max", var]))

    def get_var_min(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(" ".join(["min", var]))

    def get_var_mean(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(" ".join(["mean", var]))

    def get_var_max_slope(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_max(f"d{var}/dt")

    def get_var_min_slope(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_min(f"d{var}/dt")

    def get_var_count(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self._get_var(" ".join([var, "count"]))

    def get_event_data(
        self, name: str, type: Optional[str] = "time"
    ) -> np.ndarray[Any, np.dtype[np.float64]]:
        # event data feature names have format: "{name} event {time/var} {event idx}"
        # name - distinguish event types (up/down, localmax/localmin) in some observers, others (nhood2) just track single type
        # type - use to extract event time or event var [e.g., var = eVar value at event time]
        #
        # return: for now, force user to specify one name that makes sense, optionally type
        event_features = [
            feature_name
            for feature_name in self._feature_names
            if name in feature_name
            and "event" in feature_name
            and "count" not in feature_name
        ]
        if len(event_features) == 0:
            raise NotImplementedError(
                f"{self._observer_type} does not track {name} event {type}s!"
            )
        data = []
        for event_idx in range(0, self._op.max_event_timestamps):
            datapoint = self._get_var(f"{name} event {type} {event_idx}")
            if np.all(datapoint == 0):
                break
            data.append(datapoint)
        return np.stack(data, axis=-1).squeeze()

    # TODO deprecate function
    def get_timestamps(
        self, var: str = "event"
    ) -> np.ndarray[Any, np.dtype[np.float64]]:
        first_key = f"{var} event time 0"
        if first_key not in self._feature_names:
            raise NotImplementedError(
                f"{self._observer_type} does not track {var} event times!"
            )
        data = []
        for key_idx in range(0, self._op.max_event_timestamps):
            datapoint = self._get_var(f"{var} event time {key_idx}")
            if np.all(datapoint == 0):
                break
            datapoint = (
                datapoint[np.newaxis] if len(datapoint.shape) == 0 else datapoint
            )
            data.append(datapoint)
        return np.stack(data, axis=1).squeeze() if data else []


# TODO[mkdocs] - make consistent
class FeatureSimulator(Simulator):
    """Simulator class that stores trajectory features, computed on-the-fly

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
    t_span : tuple[float, float], optional
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

    _device_features: np.ndarray[Any, np.dtype[np.float64]] | None
    _integrator: FeatureSimulatorBase

    def __init__(
        self,
        variables: Dict[str, float],
        parameters: Dict[str, float],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        src_file: Optional[str] = None,
        rhs_equation: Optional[OpenCLRhsEquation] = None,
        supplementary_equations: Optional[List[Callable[[Any], Any]]] = None,
        stepper: Stepper = Stepper.rk4,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 10000000,
        max_store: int = 10000000,
        nout: int = 1,
        solver_parameters: Optional[SolverParams] = None,
        t_span: Tuple[float, float] = (0.0, 1000.0),
        single_precision: bool = True,
        device_type: Optional[CLDeviceType] = None,
        vendor: Optional[CLVendor] = None,
        platform_id: Optional[int] = None,
        device_id: Optional[int] = None,
        device_ids: Optional[List[int]] = None,
        observer: Observer = Observer.basic_all_variables,
        event_var: str = "",
        feature_var: str = "",
        observer_max_event_count: int = 100,  # TODO: defaults are set in two places - here and ObserverParam wrapper
        observer_max_event_timestamps: int = 0,
        observer_min_x_amp: float = 0.0,
        observer_min_imi: float = 0.0,
        observer_neighbourhood_radius: float = 0.05,
        observer_x_up_thresh: float = 0.3,
        observer_x_down_thresh: float = 0.2,
        observer_dx_up_thresh: float = 0,
        observer_dx_down_thresh: float = 0,
        observer_eps_dx: float = 0.0,
        observer_parameters: Optional[ObserverParams] = None,
    ) -> None:

        self._observer_type = observer

        event_var_idx = (
            list(variables.keys()).index(event_var) if event_var != "" else 0
        )
        feature_var_idx = (
            list(variables.keys()).index(feature_var) if feature_var != "" else 0
        )

        # can't sync yet because observer_max_event_timestamps is needed for building cl program.
        if observer_parameters is not None:
            self._op = observer_parameters
        else:
            self._op = ObserverParams(
                event_var_idx,
                feature_var_idx,
                observer_max_event_count,
                observer_max_event_timestamps,
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
            solver_parameters=solver_parameters,
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )
        # op could come after super_init if max_event_timestamps treated like trajectory max_store

    def _create_integrator(self) -> None:
        self._integrator = FeatureSimulatorBase(
            self._pi,
            self._stepper.value,
            self._observer_type.value,
            self._op,  # remove from constructor & add set_observer_pars below.
            self._single_precision,
            self._runtime,
            _clode_root_dir,
        )
        # self.set_observer_parameters()

    def set_observer(self, observer_type: Observer):
        """Change the observer"""
        if observer_type != self._observer_type:
            self._integrator.set_observer(observer_type)
            self._cl_program_is_valid = False

    # Changing solver parameters does not require re-building CL program
    def set_observer_parameters(
        self,
        op: Optional[ObserverParams] = None,
        event_var: Optional[str] = None,
        feature_var: Optional[str] = None,
        max_event_count: Optional[int] = None,
        max_event_timestamps: Optional[
            int
        ] = None,  # NOTE! Changing this invalidates the CL program!!
        min_amp: Optional[float] = None,
        min_imi: Optional[float] = None,
        nhood_radius: Optional[float] = None,
        x_up_threshold: Optional[float] = None,
        x_down_threshold: Optional[float] = None,
        dx_up_threshold: Optional[float] = None,
        dx_down_threshold: Optional[float] = None,
        eps_dx: Optional[float] = None,
    ) -> None:
        """Update any of the solver parameters and push to device

        Args:
            parameters (np.array): The parameters.

        Returns:
            None
        """
        if op is not None:
            self._op = op
        else:
            if event_var is not None:
                self._op.e_var_ix = self.variable_names.index(event_var)
            if feature_var is not None:
                self._op.f_var_ix = self.variable_names.index(feature_var)
            if max_event_count is not None:
                self._op.max_event_count = max_event_count
            if max_event_timestamps is not None:
                if self._op.max_event_timestamps != max_event_timestamps:
                    self._cl_program_is_valid = False  # NOTE! Changing max_event_timestamps invalidates the CL program!!
                self._op.max_event_timestamps = max_event_timestamps
            if min_amp is not None:
                self._op.min_amp = min_amp
            if min_imi is not None:
                self._op.min_imi = min_imi
            if nhood_radius is not None:
                self._op.nhood_radius = nhood_radius
            if x_up_threshold is not None:
                self._op.x_up_threshold = x_up_threshold
            if x_down_threshold is not None:
                self._op.x_down_threshold = x_down_threshold
            if dx_up_threshold is not None:
                self._op.dx_up_threshold = dx_up_threshold
            if dx_down_threshold is not None:
                self._op.dx_down_threshold = dx_down_threshold
            if eps_dx is not None:
                self._op.eps_dx = eps_dx
        self._integrator.set_observer_params(self._op)

    def get_observer_parameters(self):
        """Get the current observer parameter struct"""
        return self._integrator.get_observer_params()

    def get_feature_names(self) -> List[str]:
        """Get the list of feature names for the current observer"""
        return self._integrator.get_feature_names()

    def features(
        self,
        t_span: Optional[Tuple[float, float]] = None,
        initialize_observer: Optional[bool] = None,
        update_x0: bool = True,
        fetch_results: bool = True,
    ) -> Optional[ObserverOutput]:
        """Run a simulation with feature detection.

        Args:
        t_span (tuple[float, float]): Time interval for integration.
        initialize_observer (bool): Whether the observer data be initialized
        update_x0 (bool): After the simulation, whether to overwrite the initial state buffer with the final state
        fetch_results (bool): Whether to fetch the feature results from the device and return them here

        Returns:
            ObserverOutput | None
        """
        # if not self._cl_program_is_valid:
        #     self._integrator.build_cl()
        #     self._cl_program_is_valid = True

        if t_span is not None:
            self.set_tspan(t_span=t_span)

        if initialize_observer is not None:
            print(f"Setting {initialize_observer=}")
            self._integrator.features(initialize_observer)
        else:
            self._integrator.features()
        # invalidate _device_features
        # invalidate _device_final_state

        if update_x0:
            self._integrator.shift_x0()
            # invalidate _device_initial_state

        if fetch_results:
            return self.get_observer_results()

    def get_observer_results(self) -> ObserverOutput:
        """Get the features measured by the observer

        Returns:
            ObserverOutput: object containing features that summarize trajectories
        """
        self._device_features = self._integrator.get_f()
        self._num_features = self._integrator.get_n_features()

        if self._num_features is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._device_features is None:
            raise ValueError("Must run trajectory() before getting trajectory data")

        shape = (self._ensemble_size, self._num_features)
        self._device_features = np.array(
            self._device_features[: np.prod(shape)]
        ).reshape(shape, order="F")

        return ObserverOutput(
            self._op,
            self._device_features,
            self._num_features,
            self.variable_names,
            self._observer_type,
            self._integrator.get_feature_names(),
            self._ensemble_shape,
        )
