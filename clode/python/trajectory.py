import typing

from .xpp_parser import convert_xpp_file
from .stepper import Stepper
from .runtime import _get_runtime
from .runtime import _clode_root_dir
import numpy as np
from .runtime import _get_clode

_clode = _get_clode()


class CLODETrajectory:

    def __init__(
        self,
        src_file: str,
        variable_names: list[str],
        parameter_names: list[str],
        aux: typing.Optional[list[str]] = None,
        num_noise: int = 1,
        event_var: str = "",
        feature_var: str = "",
        tspan: tuple[float, float] = (0.0, 1000.0),
        stepper: Stepper = Stepper.euler,
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

        self._data = None
        self._output_trajectories = None
        self._time_steps = None
        self._number_of_simulations = None
        self._initial_conditions = None
        self._var_values = None
        self._n_stored = None
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

        self._trajectory = _clode.clode_trajectory(self._pi, stepper.value,
                                                   single_precision,
                                                   _get_runtime(),
                                                   _clode_root_dir)

        self.tspan = tspan

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

        self._data = None
        self._number_of_simulations = parameters.shape[0]

        self._trajectory.build_cl()
        self._trajectory.initialize(self.tspan,
                                    x0.transpose().flatten(),
                                    parameters.transpose().flatten(), self._sp)
        self._trajectory.seed_rng(1)

    def transient(self, update_x0=False):
        self._trajectory.transient()
        if update_x0:
            self.shift_x0()

    def shift_x0(self):
        self._trajectory.shift_x0()

    def trajectory(self):
        self._trajectory.trajectory()
        self._n_stored = self._trajectory.get_n_stored()
        self._time_steps = self._trajectory.get_t()
        self._output_trajectories = self._trajectory.get_x()
        self._initial_conditions = self._trajectory.get_x0()

    def get_time_steps(self):
        return np.array(self._time_steps[:self._n_stored[0]])

    def get_trajectory(self, simulation_id: int = 0):
        if simulation_id >= self._number_of_simulations:
            raise ValueError(
                f"Only {self._number_of_simulations} simulations were run, " +
                f"simulation id {simulation_id} not valid")
        if self._data is not None:
            return self._data

        n_stored = self._n_stored[simulation_id]

        if n_stored == 0:
            return np.array()

        shape = (self._number_of_simulations, len(self.vars), n_stored)
        arr = np.array(self._output_trajectories[:np.prod(shape)])
        self._data = arr.reshape(shape, order="F").transpose((0, 2, 1))
        return self._data
