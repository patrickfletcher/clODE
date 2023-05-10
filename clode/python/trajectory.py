from typing import List, Optional, Tuple

import numpy as np

from .runtime import _clode_root_dir, get_cpp, initialise_runtime
from .stepper import Stepper
from .xpp_parser import convert_xpp_file

_clode = get_cpp()


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
        self._number_of_simulations = None
        self._initial_conditions = None
        self._var_values = None
        self._n_stored = None
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

        self._data = None
        self._number_of_simulations = parameters.shape[0]

        if tspan is not None:
            self.tspan=tspan

        self._trajectory.initialize(
            self.tspan,
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
            self._sp,
        )
        self.seed_rng(seed)
        
    def seed_rng(self, seed: int|None = None):
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
        self._n_stored = self._trajectory.get_n_stored()
        self._time_steps = self._trajectory.get_t()
        self._output_trajectories = self._trajectory.get_x()
        self._initial_conditions = self._trajectory.get_x0()

    def get_time_steps(self):
        return np.array(self._time_steps[: self._n_stored[0]])

    def get_trajectory(self, simulation_id: int = 0):
        # TODO: simulation id not implemented?
        if simulation_id >= self._number_of_simulations:
            raise ValueError(
                f"Only {self._number_of_simulations} simulations were run, "
                + f"simulation id {simulation_id} not valid"
            )
        if self._data is not None:
            return self._data

        n_stored = self._n_stored[simulation_id]

        if n_stored == 0:
            return np.array()

        shape = (self._number_of_simulations, len(self.vars), n_stored)
        arr = np.array(self._output_trajectories[: np.prod(shape)])
        self._data = arr.reshape(shape, order="F").transpose((0, 2, 1))
        return self._data

    def print_devices(self) -> None:
        self._runtime.print_devices()

    def get_max_memory_alloc_size(self) -> int:
        return self._runtime.get_max_memory_alloc_size()

    def get_device_cl_version(self) -> str:
        return self._runtime.get_device_cl_version()

    def get_double_support(self) -> bool:
        return self._runtime.get_double_support()
