from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from .runtime import _clode_root_dir, get_cpp, initialise_runtime
from .stepper import Stepper
from .xpp_parser import convert_xpp_file

_clode = get_cpp()

# Output datastructure.  We have a collection of num_simulations trajectories, (t, X), stacked into matrices. Each trajectory may have a different number of total stored time steps.
# - for convenience, would be nice to have access patterns something like: 
# >>> trajectory_output.t[0] --> t[: nstored[0], 0]
# >>> trajectory_output.X[0] --> X[: nstored[0], :, 0]) 
# >>> trajectory_output.X[var, 0] --> X[: nstored[0], var, 0]
# >>> trajectory_output.t, .X --> (t[: max(nstored), :], X[: max(nstored), :, :]) 
class TrajectoryOutput:
    def __init__(
        self,
        number_of_simulations: int,
        output_time_steps: np.array,
        output_trajectories: np.array,
        n_stored: np.array,
        max_store: int,
        variables: list[str],
        ):

        self._number_of_simulations = number_of_simulations
        self._n_stored = n_stored
        self._vars = variables
        
        max_n_stored = max(self._n_stored)
        if max_n_stored== 0:
            return np.array()
        
        # time_steps has one column per simulation (to support adaptive steppers)
        shape = (number_of_simulations, max_store)
        arr = np.array(output_time_steps[: np.prod(shape)])
        self._time_steps = arr.reshape(shape, order="F").transpose((1, 0))
        self._time_steps = self._time_steps[: max_n_stored, :]

        shape = (number_of_simulations, len(self.vars), max_store)
        arr = np.array(output_trajectories[: np.prod(shape)])
        self._trajectory_data = arr.reshape(shape, order="F").transpose((2, 1, 0))
        self._trajectory_data = self._trajectory_data[: max_n_stored, :, :]

    # def 


class CLODETrajectory:
    _runtime: _clode.opencl_resource | None = None

    def __init__(
        self,
        src_file: str,
        variable_names: List[str],
        parameter_names: List[str],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        tspan: Tuple[float, float] = (0.0, 1000.0),
        stepper: Stepper = Stepper.rk4,
        single_precision: bool = True,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 1000000,
        max_store: int = 1000000,
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

        self._data = None
        self._output_trajectories = None
        self._time_steps = None
        self._output_time_steps = None
        self._number_of_simulations = None
        self._initial_conditions = None
        self._var_values = None
        self._n_stored = None
        self._max_store = max_store
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

        self._runtime = initialise_runtime(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )

        self._trajectory = _clode.clode_trajectory(
            self._pi, stepper.value, single_precision, self._runtime, _clode_root_dir
        )

        self.tspan = tspan
        self._trajectory.build_cl()

    def initialize(
        self,
        x0: np.array,
        parameters: np.array,
        tspan: Tuple[float, float] | None = None,
        seed: int | None = None,
    ):

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

        self._data = None
        self._time_steps = None
        self._number_of_simulations = parameters.shape[0]

        if tspan is not None:
            self.tspan = tspan

        self._trajectory.initialize(
            self.tspan,
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
            self._sp,
        )
        self.seed_rng(seed)

    def seed_rng(self, seed: int | None = None):
        if seed is not None:
            self._trajectory.seed_rng(seed)
        else:
            self._trajectory.seed_rng()

    def set_tspan(self, tspan: Tuple[float, float]):
        self.tspan = tspan
        self._trajectory.set_tspan(tspan)

    def set_problem_data(self, x0: np.array, parameters: np.array):
        self._trajectory.set_problem_data(
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
        )

    def set_x0(self, x0: np.array):
        self._trajectory.set_x0(
            x0.transpose().flatten(),
        )

    def set_parameters(self, parameters: np.array):
        self._trajectory.set_pars(
            parameters.transpose().flatten(),
        )

    def transient(self, update_x0: bool = True):
        self._trajectory.transient()
        if update_x0:
            self.shift_x0()

    def shift_tspan(self):
        self._trajectory.shift_tspan()

    def shift_x0(self):
        self._trajectory.shift_x0()

    def trajectory(self):
        self._trajectory.trajectory()

    def get_trajectory(self):
        # fetch data from device
        self._n_stored = self._trajectory.get_n_stored()
        self._output_time_steps = self._trajectory.get_t()
        self._output_trajectories = self._trajectory.get_x()
        
        # time_steps has one column per simulation (to support adaptive steppers)
        shape = (self._number_of_simulations, self._max_store)
        arr = np.array(self._output_time_steps[: np.prod(shape)])
        self._time_steps = arr.reshape(shape, order="F").transpose((1, 0))

        shape = (self._number_of_simulations, len(self.vars), self._max_store)
        arr = np.array(self._output_trajectories[: np.prod(shape)])
        self._data = arr.reshape(shape, order="F").transpose((2, 1, 0))

        # # if want to keep as matrices... need upper bound. trajectories still need to be cut to their specific n_stored...
        # max_stored = max(self._n_stored)
        # if max_stored == 0:
        #     return np.array()
        # self._time_steps = self._time_steps[: max_stored, :]
        # self._data = self._data[: max_stored, :, :]

        # alternatively, keep a list of trajectories, each stored as dict:
        result = list()
        for i in range(self._number_of_simulations):
            ni = self._n_stored[i]
            ti = self._time_steps[: ni, i]
            xi = self._data[: ni, :, i]
            result.append( {"t": ti, "X": xi} )

        return result

    def get_final_state(self):
        self._final_state = self._features.get_xf()
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
