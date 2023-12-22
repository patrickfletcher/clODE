from __future__ import annotations

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
    ) -> None:
        """Initialise a CLODE trajectory object.

        Args:
            src_file (str): The path to the source file.
            variable_names (List[str]): The names of the variables to be integrated.
            parameter_names (List[str]): The names of the parameters to be integrated.
            aux: The names of the auxiliary variables.
            num_noise: The number of noise variables.
            tspan: The time span of the integration.
            stepper: The integration method.
                * `Stepper.euler` for the Euler method.
                * `Stepper.rk4` for the 4th order Runge-Kutta method.
            single_precision: Whether to use single precision.
            dt: The initial time step.
            dtmax: The maximum time step.
            abstol: The absolute tolerance.
            reltol: The relative tolerance.
            max_steps: The maximum number of steps.
            max_store: The maximum number of steps to store.
            nout: The number of output points.
            device_type: The type of device to use.
            vendor: The vendor of the device to use.
            platform_id: The platform ID of the device to use.
            device_id: The device ID of the device to use.
            device_ids: The device IDs of the devices to use.
        """
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
    ) -> None:
        """Initialize the trajectory.

        :param x0: The initial conditions.
        :param parameters: The parameters.
        :param tspan: The time span of the integration.
        :param seed: The seed to use for the random number generator. If None,
            a random seed is used.
        """

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

    def seed_rng(self, seed: int | None = None) -> None:
        """Seed the random number generator used for noise.

        :param seed: The seed to use for the random number generator. If None,
            a random seed is used.
        """
        if seed is not None:
            self._trajectory.seed_rng(seed)
        else:
            self._trajectory.seed_rng()

    def set_tspan(self, tspan: Tuple[float, float]) -> None:
        """Set the time span of the differential equation.

        :param tspan: The time span of the differential equation.
        """
        self.tspan = tspan
        self._trajectory.set_tspan(tspan)

    def set_problem_data(self, x0: np.array, parameters: np.array) -> None:
        """Set the initial conditions and parameters of the differential equation.

        The initial conditions and parameters are flattened and transposed to
        match the CLODE convention of storing initial conditions and parameters
        in a row-major order.

        :param x0: The initial conditions of the differential equation.
        :param parameters: The parameters of the differential equation.
        """
        self._trajectory.set_problem_data(
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
        )

    def set_x0(self, x0: np.array) -> None:
        """Set the initial conditions of the differential equation.

        The initial conditions are flattened and transposed to match the CLODE
        convention of storing initial conditions in a row-major order.

        :param x0: The initial conditions of the differential equation.
        """
        self._trajectory.set_x0(
            x0.transpose().flatten(),
        )

    def set_parameters(self, parameters: np.array) -> None:
        """Set the parameters of the differential equation.

        The parameters are flattened and transposed to match the CLODE
        convention of storing parameters in a row-major order.

        :param parameters: The parameters of the differential equation.
        """
        self._trajectory.set_pars(
            parameters.transpose().flatten(),
        )

    def transient(self, update_x0: bool = True) -> None:
        """Simulate a transient and (optionally) update the initial conditions.

        :param update_x0: Whether to update the initial conditions after
            simulating a transient.
        """
        self._trajectory.transient()
        if update_x0:
            self.shift_x0()

    def shift_tspan(self) -> None:
        """Update tspan after simulating a transient."""
        self._trajectory.shift_tspan()

    def shift_x0(self) -> None:
        """Update x0 after simulating a transient."""
        self._trajectory.shift_x0()

    def trajectory(self) -> None:
        """Computes the trajectory of the differential equation."""
        self._trajectory.trajectory()
        self._n_stored = self._trajectory.get_n_stored()
        self._output_time_steps = self._trajectory.get_t()
        self._output_trajectories = self._trajectory.get_x()
        self._initial_conditions = self._trajectory.get_x0()

    def get_trajectory(self) -> (np.array, np.array, int):
        """Return the trajectory variables.

        :return:
            The time steps, the trajectory data and the number of time steps
            stored.
        """
        # For now just return all the data.
        # if simulation_id is not None and simulation_id >= self._number_of_simulations:
        #     raise ValueError(
        #         f"Only {self._number_of_simulations} simulations were run, "
        #         + f"simulation id {simulation_id} not valid"
        #     )

        if self._time_steps is not None and self._data is not None:
            return self._time_steps, self._data, self._n_stored

        max_stored = max(self._n_stored)
        if max_stored == 0:
            return np.array([]), np.array([]), 0

        # time_steps has one column per simulation (to support adaptive steppers)
        shape = (self._number_of_simulations, self._max_store)
        arr = np.array(self._output_time_steps[: np.prod(shape)])
        self._time_steps = arr.reshape(shape, order="F").transpose((1, 0))
        self._time_steps = self._time_steps[:max_stored, :]

        shape = (self._number_of_simulations, len(self.vars), self._max_store)
        arr = np.array(self._output_trajectories[: np.prod(shape)])
        self._data = arr.reshape(shape, order="F").transpose((2, 1, 0))
        self._data = self._data[:max_stored, :, :]

        return self._time_steps, self._data, self._n_stored

    def print_devices(self) -> None:
        """Print the available OpenCL devices."""
        self._runtime.print_devices()

    def get_max_memory_alloc_size(self) -> int:
        """Return the maximum memory allocation size in bytes."""
        return self._runtime.get_max_memory_alloc_size()

    def get_device_cl_version(self) -> str:
        """Return the OpenCL version of the device."""
        return self._runtime.get_device_cl_version()

    def get_double_support(self) -> bool:
        """Return whether the device supports double precision."""
        return self._runtime.get_double_support()
