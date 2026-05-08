from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from ._backends.factory import create_feature_backend
from ._backends.protocol import FeatureBackend
from .observers.types import Observer, ObserverParams
from .problem.python import OpenCLRhsEquation
from .runtime import CLDeviceType, CLVendor, _clode_root_dir
from .solver import Simulator, Stepper
from .simulation.params import SolverParams
from .simulation.results import ObserverOutput


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

    _device_features: np.ndarray[Any, np.dtype[np.float64]] | None = None
    _num_features: int | None = None
    _integrator: FeatureBackend

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
        self._integrator = create_feature_backend(
            self._pi,
            self._rhs_source,
            self._stepper.value,
            self._observer_type.value,
            self._op,
            self._single_precision,
            self._runtime,
            _clode_root_dir,
            runtime_selection=self._runtime_selection,
        )
        # self.set_observer_parameters()

    def _invalidate_feature_cache(self) -> None:
        self._device_features = None
        self._num_features = None
        self._invalidate_solution_cache()

    def set_observer(self, observer_type: Observer):
        """Change the observer"""
        if observer_type != self._observer_type:
            self._integrator.set_observer(observer_type.value)
            self._observer_type = observer_type
            self._cl_program_is_valid = False
            self._invalidate_feature_cache()

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
        """Update observer parameters and push them to the device."""
        current_max_event_timestamps = self._op.max_event_timestamps

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

        if self._op.max_event_timestamps != current_max_event_timestamps:
            self._cl_program_is_valid = False

        self._integrator.set_observer_params(self._op)
        self._invalidate_feature_cache()

    def get_observer_parameters(self):
        """Get the current observer parameter struct"""
        return self._integrator.get_observer_params()

    def get_feature_names(self) -> List[str]:
        """Get the list of feature names for the current observer"""
        return self._integrator.get_feature_names()

    def is_observer_initialized(self):
        """Get whether the current observer is initialized"""
        return self._integrator.is_observer_initialized()

    def initialize_observer(self):
        """run the observer's initialization warmup pass, if it has one"""
        self._ensure_cl_program()
        self._integrator.initialize_observer()

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
        self._ensure_cl_program()

        if t_span is not None:
            self.set_tspan(t_span=t_span)

        if initialize_observer is not None:
            self._integrator.features(initialize_observer)
        else:
            self._integrator.features()

        self._invalidate_feature_cache()

        if update_x0:
            self._integrator.shift_x0()
            # invalidate _device_initial_state
            self._device_initial_state = None

        if fetch_results:
            return self.get_observer_results()

    def get_observer_results(self) -> ObserverOutput:
        """Get the features measured by the observer

        Returns:
            ObserverOutput: object containing features that summarize trajectories
        """
        if self._device_features is None:
            self._device_features = self._integrator.get_f()
            self._num_features = self._integrator.get_n_features()

        if self._device_features is None or self._num_features is None:
            raise ValueError("Must run features() before getting observer results")

        self._device_features = np.array(
            self._device_features, dtype=np.float64
        ).reshape((self._ensemble_size, self._num_features), order="F")

        return ObserverOutput(
            self._op,
            self._device_features,
            self._num_features,
            self.variable_names,
            self._observer_type,
            self._integrator.get_feature_names(),
            self._ensemble_shape,
        )
