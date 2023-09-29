from typing import List, Optional, Tuple

import numpy as np

from .observer import Observer, ObserverOutput
from .runtime import _clode_root_dir, get_cpp, initialise_runtime
from .stepper import Stepper
from .xpp_parser import convert_xpp_file

_clode = get_cpp()


class CLODEFeatures:
    """
    A class for computing features from a CLODE model.

    Parameters
    ----------
    src_file : str
        Path to the CLODE model source file.
    variable_names : list[str]
        List of variable names in the model.
    parameter_names : list[str]
        List of parameter names in the model.
    aux : list[str], optional
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
    >>> model = clode.CLODEFeatures(
    ...     src_file="examples/lorenz96.c",
    ...     variable_names=["x"],
    ...     parameter_names=["F"],

    ... )
    >>> model.set_parameter_values({"F": 8.0})
    >>> model.set_initial_values({"x": np.random.rand(40)})
    >>> model.simulate()
    >>> model.plot()
    >>> plt.show()"""

    _runtime: _clode.opencl_resource | None = None

    def __init__(
        self,
        src_file: str,
        variable_names: List[str],
        parameter_names: List[str],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        event_var: str = "",
        feature_var: str = "",
        observer_max_event_count: int = 100,
        observer_min_x_amp: float = 1.0,
        observer_min_imi: float = 1,
        observer_neighbourhood_radius: float = 0.01,
        observer_x_up_thresh: float = 0.3,
        observer_x_down_thresh: float = 0.2,
        observer_dx_up_thresh: float = 0,
        observer_dx_down_thresh: float = 0,
        observer_eps_dx: float = 1e-7,
        tspan: Tuple[float, float] = (0.0, 1000.0),
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
        device_type: _clode.cl_device_type | None = None,
        vendor: _clode.cl_vendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
        device_ids: List[int] | None = None,
    ):
        if src_file.endswith(".xpp"):
            input_file = convert_xpp_file(src_file)
        else:
            input_file = src_file

        self._final_state = None
        self._num_result_features = None
        self._result_features = None
        if aux is None:
            aux = []

        self.vars = variable_names
        self.pars = parameter_names
        self.aux_variables = aux
        self._pi = _clode.problem_info(
            input_file,
            len(variable_names),
            len(parameter_names),
            len(aux),
            num_noise,
            variable_names,
            parameter_names,
            aux,
        )
        self._sp = _clode.solver_params(
            dt, dtmax, abstol, reltol, max_steps, max_store, nout
        )

        event_var_idx = variable_names.index(event_var) if event_var != "" else 0
        feature_var_idx = variable_names.index(feature_var) if feature_var != "" else 0

        self._op = _clode.observer_params(
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

        self._runtime = initialise_runtime(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )

        self._features = _clode.clode_features(
            self._pi,
            stepper.value,
            observer.value,
            single_precision,
            self._runtime,
            _clode_root_dir,
        )

        self.tspan = tspan
        self._observer_type = observer
        self._features.build_cl()

    def initialize(self, x0: np.array, parameters: np.array, tspan:Tuple[float, float]|None = None, seed: int|None = None):

        if len(x0.shape) != 2:
            raise ValueError("Must provide rows of initial variables")

        if x0.shape[1] != len(self.vars):
            raise ValueError(
                f"Length of initial condition vector {len(x0.shape[1])}"
                f" does not match number of variables {len(self.vars)}"
            )

        if len(parameters.shape) != 2:
            raise ValueError("Must provide rows of parameters")

        if parameters.shape[1] != len(self.pars):
            raise ValueError(
                f"Length of parameters vector {parameters.shape[1]}"
                f" does not match number of parameters {len(self.pars)}"
            )
        
        if tspan is not None:
            self.tspan=tspan

        self._features.initialize(
            self.tspan,
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
            self._sp,
            self._op,
        )
        self.seed_rng(seed)

    def seed_rng(self, seed: int|None = None):
        if seed is not None:
            self._features.seed_rng(seed)
        else:
            self._features.seed_rng()
        
    def set_tspan(self, new_tspan: Tuple[float, float]):
        self.tspan = new_tspan
        self._features.set_tspan(new_tspan)

    def set_problem_data(self, x0: np.array, parameters: np.array):
        self._features.set_problem_data(
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
            )
        
    def set_x0(self, x0: np.array):
        self._features.set_x0(
            x0.transpose().flatten(),
            )
        
    def set_parameters(self, parameters: np.array):
        self._features.set_pars(
            parameters.transpose().flatten(),
            )

    def transient(self, update_x0: bool = True):
        self._features.transient()
        if update_x0:
            self.shift_x0()

    def shift_tspan(self):
        self._features.shift_tspan()

    def shift_x0(self):
        self._features.shift_x0()

    def features(self, initialize_observer: Optional[bool] = None):
        if initialize_observer is not None:
            print("Reinitializing observer")
            self._features.features(initialize_observer)
        else:
            self._features.features()
        self._result_features = self._features.get_f()
        self._num_result_features = self._features.get_n_features()
        self._final_state = self._features.get_xf()

    def get_observer_results(self):
        return ObserverOutput(
            self._op,
            np.array(self._result_features),
            self._num_result_features,
            self.vars,
            self._observer_type,
            self._features.get_feature_names(),
        )

    def get_final_state(self):
        final_state = np.array(self._final_state)
        return final_state.reshape(
            (len(self.vars), len(final_state) // len(self.vars))
        ).transpose()

    def print_devices(self) -> None:
        self._runtime.print_devices()

    def get_max_memory_alloc_size(self) -> int:
        return self._runtime.get_max_memory_alloc_size()

    def get_device_cl_version(self) -> str:
        return self._runtime.get_device_cl_version()

    def get_double_support(self) -> bool:
        return self._runtime.get_double_support()
