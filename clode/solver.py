from __future__ import annotations

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from .function_converter import OpenCLConverter, OpenCLRhsEquation
from .runtime import (
    CLDeviceType,
    CLVendor,
    OpenCLResource,
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

    _integrator: SimulatorBase
    _pi: ProblemInfo
    _sp: SolverParams
    _runtime: OpenCLResource
    t_span: Tuple[float, float]
    _variables: Optional[Dict[str, np.ndarray]] = None
    _variable_defaults: Dict[str, float]
    _parameters: Optional[Dict[str, np.ndarray]] = None
    _parameter_defaults: Dict[str, float]
    _ensemble_size: int
    _single_precision: bool
    _stepper: Stepper

    @property
    def variable_names(self) -> List[str]:
        """Get the variable names.

        Returns:
            List[str]: The variable names.
        """
        return list(self._variable_defaults.keys())

    @property
    def parameter_names(self) -> List[str]:
        """Get the parameter names.

        Returns:
            List[str]: The parameter names.
        """
        return list(self._parameter_defaults.keys())

    @property
    def is_initialized(self) -> bool:
        """Get whether the simulator is initialized.

        Returns:
            bool: Whether the simulator is initialized.
        """
        # return self._integrator.is_initialized()
        return True

    def _find_ensemble_size(
        self,
        variables: Optional[
            Dict[str, Union[float, np.ndarray[np.dtype[np.float64]], List[float]]]
        ],
        parameters: Optional[
            Dict[str, Union[float, np.ndarray[np.dtype[np.float64]], List[float]]]
        ],
    ) -> int:
        """Find the length of the arrays to be passed to the OpenCL kernel.

        Args:
            variables (Dict[str, Union[float, np.ndarray, List[float]]]): The variables.
            parameters (Dict[str, Union[float, np.ndarray, List[float]]]): The parameters.

        Returns:
            int: The length of the arrays to be passed to the OpenCL kernel.
        """
        array_length: Optional[int] = None

        if variables is None:
            array_length = 1
        else:
            for key, value in variables.items():
                if isinstance(value, (np.ndarray, List)):
                    if array_length is None or array_length == 1:
                        array_length = len(value)
                    elif array_length != len(value):
                        raise ValueError(
                            f"Variable {key} has length {len(value)} "
                            f"but previous variables have length {array_length}"
                        )
                elif isinstance(value, float):
                    if array_length is None:
                        array_length = 1

        if parameters is None:
            array_length = 1 if array_length is None else array_length
        else:
            for key, value in parameters.items():
                if isinstance(value, (np.ndarray, List)):
                    if array_length is None or array_length == 1:
                        array_length = len(value)
                    elif array_length != len(value):
                        raise ValueError(
                            f"Parameter {key} has length {len(value)} "
                            f"but previous parameters have length {array_length}"
                        )
                elif isinstance(value, float):
                    if array_length is None:
                        array_length = 1
        if array_length is None:
            raise ValueError("No variables or parameters specified")
        return array_length

    def _create_cl_arrays(
        self, data: Dict[str, Union[float, np.ndarray, List[float]]], array_length: int
    ) -> Dict[str, np.ndarray]:
        cl_data: Dict[str, np.ndarray] = {}
        for key, value in data.items():
            array: np.ndarray

            if isinstance(value, float):
                array = np.full(array_length, value)
            elif isinstance(value, np.ndarray):
                array = value
            elif isinstance(value, List):
                array = np.array(value)
            else:
                raise ValueError(f"Invalid type for {key}")
            # Check that array is the correct length
            if len(array) != array_length:
                raise ValueError(
                    f"Array {key} has length {len(array)} but should have length {array_length}"
                )
            cl_data[key] = array
        return cl_data

    def _handle_clode_rhs_cl_file(
        self,
        src_file: str | None = None,
        rhs_equation: OpenCLRhsEquation | None = None,
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
            eqn = converter.convert_to_opencl(
                rhs_equation, mutable_args=[3, 4], function_name="getRHS"
            )
            input_file = "clode_rhs.cl"
            with open(input_file, "w") as ff:
                ff.write(eqn)
        else:
            raise ValueError("Must specify either src_file or rhs_equation")

        return input_file

    def __init__(
        self,
        variables: Dict[str, float],
        parameters: Dict[str, float],
        src_file: str | None = None,
        rhs_equation: OpenCLRhsEquation | None = None,
        supplementary_equations: List[Callable[[Any], Any]] | None = None,
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        t_span: Tuple[float, float] = (
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
        self._stepper = stepper
        self._single_precision = single_precision

        self._variable_defaults = variables
        self._parameter_defaults = parameters
        self._ensemble_size = 1

        input_file = self._handle_clode_rhs_cl_file(
            src_file, rhs_equation, supplementary_equations
        )

        if aux is None:
            aux = []

        self._max_store = max_store

        self.aux_variables = aux
        self._pi = ProblemInfo(
            input_file,
            self.variable_names,
            self.parameter_names,
            aux,
            num_noise,
        )
        self._sp = SolverParams(dt, dtmax, abstol, reltol, max_steps, max_store, nout)

        self.tspan = t_span

        # _runtime as an instance variable
        self._runtime = initialize_runtime(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )

        self._build_integrator()

        self._integrator.build_cl()

        self._init_integrator()

    def _pack_data(self) -> Tuple[np.ndarray, np.ndarray]:
        """Pack the data into a tuple.

        Returns:
            Tuple[np.ndarray, np.ndarray]: The packed data.
        """
        # Pack the varibles, transpose and flatten
        # so that they are in the correct format for the OpenCL kernel
        if self._variables is not None:
            vars_array = np.array(list(self._variables.values()))
        else:
            vars_array = np.array(list(self._variable_defaults.values()))

        if self._parameters is not None:
            pars_array = np.array(list(self._parameters.values()))
        else:
            pars_array = np.array(list(self._parameter_defaults.values()))

        return vars_array.flatten(), pars_array.flatten()

    def set_repeat_ensemble(self, num_repeats: int) -> None:
        """Set the number of repeats for the ensemble.

        Args:
            num_repeats (int): The number of repeats for the ensemble.

        Returns:
            None
        """
        self._ensemble_size = num_repeats

        self._variables = self._create_cl_arrays(self._variable_defaults, num_repeats)

        self._parameters = self._create_cl_arrays(self._parameter_defaults, num_repeats)

        # if self.is_initialized:
        #     vars_array, pars_array = self._pack_data()
        #     self._integrator.set_problem_data(vars_array, pars_array)
        # else:
        self._init_integrator()

    def set_ensemble(
        self,
        variables: Optional[
            Union[
                np.ndarray[np.dtype[np.float64]],
                Dict[str, Union[np.ndarray[np.dtype[np.float64]], List[float]]],
            ]
        ] = None,
        parameters: Optional[
            Union[
                np.ndarray[np.dtype[np.float64]],
                Dict[str, Union[np.ndarray[np.dtype[np.float64]], List[float]]],
            ]
        ] = None,
        seed: Optional[int] = None,
    ) -> None:
        if isinstance(variables, np.ndarray):
            # We test that the array is the correct length
            if len(variables.shape) != 2:
                raise ValueError("Variables must be a matrix")
            if variables.shape[1] != len(self.variable_names):
                raise ValueError(
                    f"Variables must have {len(self.variable_names)} columns"
                )
            variables = {
                key: variables[:, index]
                for index, key in enumerate(self.variable_names)
            }
        if isinstance(parameters, np.ndarray):
            # We test that the array is the correct length
            if len(parameters.shape) != 2:
                raise ValueError("Parameters must be a matrix")
            if parameters.shape[1] != len(self.parameter_names):
                raise ValueError(
                    f"Parameters must have {len(self.parameter_names)} columns"
                )
            parameters = {
                key: parameters[:, index]
                for index, key in enumerate(self.parameter_names)
            }

        cl_array_length = self._find_ensemble_size(variables, parameters)

        # Implicitly create arrays of the correct length
        # Discard self._variables and self._parameters
        local_variables: Dict[
            str, Union[float, np.ndarray[np.dtype[np.float64]], List[float]]
        ] = {}

        # Keys can be missing in variables but not in self._variable_defaults
        if variables is not None:
            for key in variables.keys():
                if key not in self._variable_defaults:
                    raise KeyError(f"Key {key} not in ODE variable defaults!")

        for key in self._variable_defaults.keys():
            if variables is not None and key in variables:
                local_variables[key] = variables[key]
            else:
                local_variables[key] = self._variable_defaults[key]

        local_parameters: Dict[
            str, Union[float, np.ndarray[np.dtype[np.float64]], List[float]]
        ] = {}

        # Keys can be missing in parameters but not in self._parameter_defaults
        if parameters is not None:
            for key in parameters.keys():
                if key not in self._parameter_defaults:
                    raise KeyError(f"Key {key} not in ODE parameter defaults!")

        for key in self._parameter_defaults.keys():
            if parameters is not None and key in parameters:
                local_parameters[key] = parameters[key]
            else:
                local_parameters[key] = self._parameter_defaults[key]

        self._variables = self._create_cl_arrays(local_variables, cl_array_length)

        self._parameters = self._create_cl_arrays(local_parameters, cl_array_length)

        self._ensemble_size = cl_array_length

        # if self.is_initialized:
        #     vars_array, pars_array = self._pack_data()
        #     self._integrator.set_problem_data(vars_array, pars_array)
        # else:
        self._init_integrator()

        # self.seed_rng(seed)

    def _build_integrator(self) -> None:
        self._integrator = SimulatorBase(
            self._pi,
            self._stepper.value,
            self._single_precision,
            self._runtime,
            _clode_root_dir,
        )

    def _init_integrator(self) -> None:
        # FeaturesSimulator overrides this
        vars_array, pars_array = self._pack_data()
        self._integrator.initialize(
            self.tspan,
            vars_array,
            pars_array,
            self._sp,
        )

    def set_initial_conditions(
        self, initial_conditions: dict[str, Union[float, np.ndarray, List[float]]]
    ) -> None:
        """Set the initial conditions.

        Args:
            initial_conditions (dict[str, Union[float, np.ndarray, List[float]]]): The initial conditions.

        Returns:
            None
        """
        self._initial_conditions = initial_conditions

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
        # self._integrator.set_problem_data(
        #     x0.transpose().flatten(),
        #     parameters.transpose().flatten(),
        # )
        raise NotImplementedError

    def set_x0(self, x0: np.array) -> None:
        """Set the initial conditions.

        Args:
            x0 (np.array): The initial conditions.

        Returns:
            None
        """
        # self._integrator.set_x0(
        #     x0.transpose().flatten(),
        # )
        raise NotImplementedError

    def set_parameters(self, parameters: np.array) -> None:
        """Set the parameters.

        Args:
            parameters (np.array): The parameters.

        Returns:
            None
        """
        # self._integrator.set_pars(
        #     parameters.transpose().flatten(),
        # )
        raise NotImplementedError

    def set_solver_parameters(
        self,
    ) -> None:
        """Set the solver parameters.

        Args:
            parameters (np.array): The parameters.

        Returns:
            None
        """
        raise NotImplementedError

    def transient(self, update_x0: bool = True) -> None:
        """Run a transient simulation.

        Args:
            update_x0 (bool, optional): Whether to update the initial conditions. Defaults to True.

        Returns:
            None
        """
        if not self.is_initialized:
            raise RuntimeError("Simulator is not initialized")
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
            (len(self._variables), len(initial_state) // len(self._variables))
        ).transpose()

    def get_final_state(self):
        """Get the final state of the simulation.

        Returns:
            np.array: The final state of the simulation.
        """
        self._final_state = self._integrator.get_xf()
        final_state = np.array(self._final_state)
        return final_state.reshape(
            (len(self.variable_names), len(final_state) // len(self.variable_names))
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
