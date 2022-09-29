import src.clode as _clode
from enum import Enum
import os
import typing

import numpy as np

_runtime = None
_clode_root_dir: str = os.path.join(os.path.dirname(os.path.dirname(__file__)), "src", "")


def _get_runtime():
    global _runtime
    if _runtime is None:
        _runtime = _clode.opencl_resource()
    return _runtime


class Stepper(Enum):
    euler = "euler"


class Observer(Enum):
    basic = "basic"


class ProblemInfo:

    def __init__(self,
                 src_file: str,
                 vars: [str],
                 pars: [str],
                 aux: [str],
                 num_noise: int):
        _pi = _clode.problem_info(src_file,
                                  len(vars),
                                  len(pars),
                                  len(aux),
                                  num_noise,
                                  vars,
                                  pars,
                                  aux)


class SolverParams:
    def __init__(self,
                 dt: float = 0.1,
                 dtmax: float = 1.0,
                 abstol: float = 1e-6,
                 reltol: float = 1e-3,
                 max_steps: int = 10000000,
                 max_store: int = 10000000,
                 nout: int = 50):
        self._sp = _clode.solver_params(dt, dtmax, abstol, reltol,
                                        max_steps, max_store, nout)


class CLODEFeatures:

    def __init__(self,
                 src_file: str,
                 variables: list[str],
                 parameters: list[str],
                 aux: typing.Optional[list[str]] = None,
                 num_noise: int = 1,
                 event_var: str = "",
                 feature_var: str = "",
                 observer_max_event_count: int = 100,
                 observer_min_x_amp: float = 1.0,
                 observer_min_imi: float = 1.0,
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
        self._final_state = None
        self._num_result_features = None
        self._result_features = None
        if aux is None:
            aux = []

        self.vars = variables
        self.pars = parameters
        self.aux_variables = aux
        self._pi = _clode.problem_info(src_file,
                                       len(variables),
                                       len(parameters),
                                       len(aux),
                                       num_noise,
                                       variables,
                                       parameters,
                                       aux)
        self._sp = _clode.solver_params(dt, dtmax, abstol, reltol,
                                        max_steps, max_store, nout)

        event_var_idx = variables.index(event_var) if event_var != "" else 0
        feature_var_idx = variables.index(feature_var) if feature_var != "" else 0

        self._op = _clode.observer_params(event_var_idx,
                                          feature_var_idx,
                                          observer_max_event_count,
                                          observer_min_x_amp,
                                          observer_min_imi,
                                          observer_neighbourhood_radius,
                                          observer_x_up_thresh,
                                          observer_x_down_thresh,
                                          observer_dx_up_thresh,
                                          observer_dx_down_thresh,
                                          observer_eps_dx)

        self._features = _clode.clode_features(self._pi,
                                               stepper.value,
                                               observer.value,
                                               single_precision,
                                               _get_runtime(),
                                               _clode_root_dir)

        self.tspan = tspan

    def initialize(self,
                   x0: np.array,
                   parameters: np.array):

        if len(x0.shape) != 2:
            raise ValueError("Must provide rows of initial variables")

        if x0.shape[1] != len(self.vars):
            raise ValueError(f"Length of initial condition vector {len(x0.shape[1])}"
                             f" does not match number of variables {len(self.vars)}")

        if len(parameters.shape) != 2:
            raise ValueError("Most provide rows of parameters")

        if parameters.shape[1] != len(self.pars):
            raise ValueError(f"Length of parameters vector {parameters.shape[1]}"
                             f" does not match number of parameters {len(self.pars)}")

        print("Init with x0", x0.transpose().flatten().shape)
        print("Init with pars", parameters.transpose().flatten().shape)
        self._features.initialize(self.tspan,
                                  x0.transpose().flatten(),
                                  parameters.transpose().flatten(),
                                  self._sp,
                                  self._op)

    def transient(self):
        self._features.transient()

    def features(self):
        self._features.features()
        self._result_features = self._features.get_f()
        self._num_result_features = self._features.get_n_features()
        self._final_state = self._features.getXf()

    def get_observer_results(self):
        # if self._result_features is None:
        #     results = self._features.features()
        #
        # else:
        #     results = self._result_features
        print("Len results", len(self._result_features), self._num_result_features)
        results = np.array(self._result_features)
        print(len(results), results.shape[0])
        print((self._num_result_features, results.shape[0] // len(self.vars)))

        return results.reshape((self._num_result_features, len(results) // self._num_result_features)).transpose()

    def get_final_state(self):
        final_state = np.array(self._final_state)
        return final_state.reshape((len(self.vars), len(final_state) // len(self.vars))).transpose()
