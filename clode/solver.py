from __future__ import annotations

from enum import Enum
from typing import Any, Callable, List, Optional, Tuple

import numpy as np

from .function_converter import OpenCLConverter, OpenCLRhsEquation
from .runtime import (
    CLDeviceType,
    CLVendor,
    ProblemInfo,
    SimulatorBase,
    SolverParams,
    _clode_root_dir,
    initialize_runtime,
)
from .xpp_parser import convert_xpp_file


class Stepper(Enum):
    euler = "euler"
    heun = "heun"
    rk4 = "rk4"
    bs23 = "bs23"
    dormand_prince = "dopri5"
    stochastic_euler = "seuler"


# base solver class with only transient()
class Simulator:
    """Base class for all simulators.

    This class is used to simulate transient behavior of a system of ODEs.

    It can be used directly to find separatrix solutions
    and variable steady states,
    or it can be used as a base class for other simulators."""

    # cleaner interface? Use the problem_info and solver_params "structs"

    def _handle_clode_rhs_cl_file(
        self,
        src_file: str | None = None,
        rhs_equation: OpenCLRhsEquation | None = None,
        mutable_args: List[str] | None = None,
        supplementary_equations: List[Callable[[Any], Any]] | None = None,
    ) -> str:
        input_file: str

        if src_file is not None and rhs_equation is not None:
            raise ValueError("Cannot specify both src_file and rhs_equation")
        elif src_file is not None:
            if src_file.endswith(".xpp"):
                input_file = convert_xpp_file(src_file)
            else:
                input_file = src_file
        elif rhs_equation is not None:
            # Convert the rhs_equation to a string
            # and write it to a file
            # using function_converter
            converter = OpenCLConverter()
            if supplementary_equations is not None:
                for eq in supplementary_equations:
                    converter.convert_to_opencl(eq)
            # TODO propagage mutable_args from __init__ to here
            if mutable_args is None:
                mutable_args = ["derivatives", "aux"]
            eqn = converter.convert_to_opencl(rhs_equation, mutable_args=mutable_args)
            input_file = "clode_rhs.cl"
            with open(input_file, "w") as ff:
                ff.write(eqn)
        else:
            raise ValueError("Must specify either src_file or rhs_equation")

        return input_file

    def __init__(
        self,
        variable_names: List[str],
        parameter_names: List[str],
        src_file: str | None = None,
        rhs_equation: OpenCLRhsEquation | None = None,
        mutable_args: List[str] | None = None,
        supplementary_equations: List[Callable[[Any], Any]] | None = None,
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        tspan: Tuple[float, float] = (
            0.0,
            1000.0,
        ),  # tspan <- realistically usually would set this as arg during integrate
        stepper: Stepper = Stepper.rk4,  # stepper <- goes with solver_params?
        single_precision: bool = True,  # precision
        dt: float = 0.1,  # solver_params
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 1000000,
        max_store: int = 1000000,
        nout: int = 1,
        device_type: CLDeviceType | None = None,  # device selection
        vendor: CLVendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
        device_ids: List[int] | None = None,
    ) -> None:
        input_file = self._handle_clode_rhs_cl_file(
            src_file, rhs_equation, mutable_args, supplementary_equations
        )

        # self._data = None
        # self._output_trajectories = None
        # self._time_steps = None
        # self._output_time_steps = None
        # self._number_of_simulations = None
        # self._initial_conditions = None
        # self._var_values = None
        # self._n_stored = None
        self._max_store = max_store
        if aux is None:
            aux = []

        self.vars = variable_names
        self.pars = parameter_names
        self.aux_variables = aux
        self._pi = ProblemInfo(
            input_file,
            variable_names,
            parameter_names,
            aux,
            num_noise,
        )
        self._sp = SolverParams(dt, dtmax, abstol, reltol, max_steps, max_store, nout)

        self.tspan = tspan

        # _runtime as an instance variable
        self._runtime = initialize_runtime(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )

        # Check whether this is being called by a derived class:
        if type(self) is Simulator:
            self._integrator = SimulatorBase(
                self._pi,
                stepper.value,
                single_precision,
                self._runtime,
                _clode_root_dir,
            )
            self._integrator.build_cl()

    # same for clODE, clODETrajectory, overridden by clODEFeatures....
    def initialize(
        self,
        x0: np.array,
        parameters: np.array,
        tspan: Tuple[float, float] | None = None,
        seed: int | None = None,
    ) -> None:
        """Initialize the clode object.

        Args:
            x0 (np.array): The initial conditions.
            parameters (np.array): The parameters.
            tspan (Tuple[float, float], optional): The time span to simulate over. Defaults to None.
            seed (int, optional): The seed for the random number generator. Defaults to None.

        Raises:
            ValueError: If the initial conditions or parameters are not the correct shape.

        Returns:
            None
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

        self._number_of_simulations = parameters.shape[0]

        if tspan is not None:
            self.tspan = tspan

        self._integrator.initialize(
            self.tspan,
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
            self._sp,
        )
        self.seed_rng(seed)

    def seed_rng(self, seed: int | None = None) -> None:
        """Seed the random number generator.

        Args:
            seed (int, optional): The seed for the random number generator. Defaults to None.

        Returns:
            None
        """

        if seed is not None:
            self._integrator.seed_rng(seed)
        else:
            self._integrator.seed_rng()

    def set_tspan(self, tspan: Tuple[float, float]) -> None:
        """Set the time span to simulate over.

        Args:
            tspan (Tuple[float, float]): The time span to simulate over.

        Returns:
            None
        """
        self.tspan = tspan
        self._integrator.set_tspan(tspan)

    def set_problem_data(self, x0: np.array, parameters: np.array) -> None:
        """Set the problem data.

        Args:
            x0 (np.array): The initial conditions.
            parameters (np.array): The parameters.

        Returns:
            None
        """
        self._integrator.set_problem_data(
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
        )

    def set_x0(self, x0: np.array) -> None:
        """Set the initial conditions.

        Args:
            x0 (np.array): The initial conditions.

        Returns:
            None
        """
        self._integrator.set_x0(
            x0.transpose().flatten(),
        )

    def set_parameters(self, parameters: np.array) -> None:
        """Set the parameters.

        Args:
            parameters (np.array): The parameters.

        Returns:
            None
        """
        self._integrator.set_pars(
            parameters.transpose().flatten(),
        )

    def set_solver_parameters(
        self,
    ) -> None:
        """Set the solver parameters.

        Args:
            parameters (np.array): The parameters.

        Returns:
            None
        """
        pass

    def transient(self, update_x0: bool = True) -> None:
        """Run a transient simulation.

        Args:
            update_x0 (bool, optional): Whether to update the initial conditions. Defaults to True.

        Returns:
            None
        """
        self._integrator.transient()
        if update_x0:
            self.shift_x0()

    def shift_tspan(self) -> None:
        """Shift the time span to the current time plus the time period.

        Returns:
            None
        """
        self._integrator.shift_tspan()

    def shift_x0(self) -> None:
        """Shift the initial conditions to the current variable values.

        Returns:
            None
        """
        self._integrator.shift_x0()

    def get_initial_state(self):
        """Get the final state of the simulation.

        Returns:
            np.array: The final state of the simulation.
        """
        self._initial_state = self._integrator.get_x0()
        initial_state = np.array(self._initial_state)
        return initial_state.reshape(
            (len(self.vars), len(initial_state) // len(self.vars))
        ).transpose()

    def get_final_state(self):
        """Get the final state of the simulation.

        Returns:
            np.array: The final state of the simulation.
        """
        self._final_state = self._integrator.get_xf()
        final_state = np.array(self._final_state)
        return final_state.reshape(
            (len(self.vars), len(final_state) // len(self.vars))
        ).transpose()

    def get_available_steppers(self) -> List[str]:
        """Get the available time steppers.

        Returns:
            List[str]
        """
        return self._integrator.get_available_steppers()
    
    def get_program_string(self) -> str:
        """Get the program string.

        Returns:
            str
        """
        return self._integrator.get_program_string()

    def print_status(self) -> None:
        """Print the simulator status info.

        Returns:
            None
        """
        self._integrator.print_status()

    def print_devices(self) -> None:
        """Print the available devices.

        Returns:
            None
        """
        self._runtime.print_devices()

    def get_max_memory_alloc_size(self) -> int:
        """Get the maximum memory allocation size.

        Returns:
            int: The maximum memory allocation size.
        """
        return self._runtime.get_max_memory_alloc_size()

    def get_device_cl_version(self) -> str:
        """Get the device OpenCL version.

        Returns:
            str: The device OpenCL version.
        """
        return self._runtime.get_device_cl_version()

    def get_double_support(self) -> bool:
        """Get whether double precision is supported.

        Returns:
            bool: Whether double precision is supported.
        """
        return self._runtime.get_double_support()
