import typing

from .xpp_parser import convert_xpp_file
from .stepper import Stepper
from .observer import Observer, ObserverOutput
from .runtime import _get_clode, _get_runtime, _clode_root_dir
import numpy as np

_clode = _get_clode()


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
    >>> plt.show()
"""
    def __init__(
        self,
        src_file: str,
        variable_names: list[str],
        parameter_names: list[str],
        aux: typing.Optional[list[str]] = None,
        num_noise: int = 1,
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
        tspan: tuple[float, float] = (0.0, 1000.0),
        stepper: Stepper = Stepper.euler,
        observer: Observer = Observer.basic,
        single_precision: bool = False,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 10000000,
        max_store: int = 10000000,
        nout: int = 50,
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
        self._pi = _clode.problem_info(input_file, len(variable_names),
                                       len(parameter_names), len(aux),
                                       num_noise, variable_names,
                                       parameter_names, aux)
        self._sp = _clode.solver_params(dt, dtmax, abstol, reltol, max_steps,
                                        max_store, nout)

        event_var_idx = variable_names.index(
            event_var) if event_var != "" else 0
        feature_var_idx = variable_names.index(
            feature_var) if feature_var != "" else 0

        self._op = _clode.observer_params(
            event_var_idx, feature_var_idx, observer_max_event_count,
            observer_min_x_amp, observer_min_imi,
            observer_neighbourhood_radius, observer_x_up_thresh,
            observer_x_down_thresh, observer_dx_up_thresh,
            observer_dx_down_thresh, observer_eps_dx)

        self._features = _clode.clode_features(self._pi, stepper.value,
                                               observer.value,
                                               single_precision,
                                               _get_runtime(), _clode_root_dir)

        self.tspan = tspan
        self._observer_type = observer

    def initialize(self, x0: np.array, parameters: np.array):

        if len(x0.shape) != 2:
            raise ValueError("Must provide rows of initial variables")

        if x0.shape[1] != len(self.vars):
            raise ValueError(
                f"Length of initial condition vector {len(x0.shape[1])}"
                f" does not match number of variables {len(self.vars)}")

        if len(parameters.shape) != 2:
            raise ValueError("Most provide rows of parameters")

        if parameters.shape[1] != len(self.pars):
            raise ValueError(
                f"Length of parameters vector {parameters.shape[1]}"
                f" does not match number of parameters {len(self.pars)}")

        self._features.build_cl()
        self._features.initialize(self.tspan,
                                  x0.transpose().flatten(),
                                  parameters.transpose().flatten(), self._sp,
                                  self._op)
        self._features.seed_rng(1)

    def transient(self, update_x0=True):
        self._features.transient()
        if update_x0:
            self.shift_x0()

    def shift_x0(self):
        self._features.shift_x0()

    def features(self, initialize_observer: typing.Optional[bool] = None):
        if initialize_observer is not None:
            print("Reinitializing observer")
            self._features.features(initialize_observer)
        else:
            self._features.features()
        self._result_features = self._features.get_f()
        self._num_result_features = self._features.get_n_features()
        self._final_state = self._features.getXf()

    def get_observer_results(self):
        return ObserverOutput(self._op, np.array(self._result_features),
                              self._num_result_features, self.vars,
                              self._observer_type,
                              self._features.get_feature_names())

    def get_final_state(self):
        final_state = np.array(self._final_state)
        return final_state.reshape(
            (len(self.vars), len(final_state) // len(self.vars))).transpose()
