from __future__ import annotations

from enum import Enum
from typing import Any, Dict, List, Literal, Callable, Optional, Tuple, Union, Mapping

import numpy as np

from .function_converter import OpenCLConverter, OpenCLRhsEquation
from .runtime import (
    CLDeviceType,
    CLVendor,
    OpenCLResource,
    _clode_root_dir,
    initialize_runtime,
    get_log_level,
    set_log_level,
    LogLevel,
)
from clode.cpp.clode_cpp_wrapper import SolverParams, ProblemInfo, SimulatorBase

from .xpp_parser import convert_xpp_file


# TODO[API]: different steppers use different parameters subsets.
# Model with classes/mixins?
# - StepperBase + mixins --> each stepper as python class.
class Stepper(Enum):
    euler = "euler"
    heun = "heun"
    rk4 = "rk4"
    bs23 = "bs23"
    dormand_prince = "dopri5"
    stochastic_euler = "seuler"


# TODO[API]: defaults for solverParams are copied in each constructor plus wrapper.
# Should be only ONE place globally.
# - Prefer the struct defaults?


class Simulator:
    """Base class for simulating an ensemble of instances of an ODE system.

    It provides the core functionality for advancing the simulation in time without
    storing any intermediate state. May be used directly when only the final state is
    of interest, or as a base class for other simulators.
    """

    _integrator: SimulatorBase
    _runtime: OpenCLResource
    _single_precision: bool
    _stepper: Stepper
    _pi: ProblemInfo

    # changes to the above items require rebuilding the CL program.
    # Flag them and rebuild if necessary on simulation function call
    _cl_program_is_valid: bool = False

    _sp: SolverParams
    _t_span: Tuple[float, float]

    _variable_defaults: Dict[str, float]
    _parameter_defaults: Dict[str, float]

    _variables: Optional[Dict[str, np.ndarray]] = None
    _parameters: Optional[Dict[str, np.ndarray]] = None

    _ensemble_size: int  # C++ layer: nPts
    _ensemble_shape: Tuple

    # 2D array shape (ensemble_size, num_parameters)
    _device_parameters: Optional[np.ndarray] = None

    # 2D array shape (ensemble_size, num_variables)
    _device_initial_state: Optional[np.ndarray] = None
    _device_final_state: Optional[np.ndarray] = None

    @property
    def variable_names(self) -> List[str]:
        """The list of ODE variable names"""
        return self._pi.vars

    @property
    def num_variables(self) -> int:
        """The number of ODE state variables"""
        return self._pi.num_var

    @property
    def parameter_names(self) -> List[str]:
        """The list of ODE system parameter names"""
        return self._pi.pars

    @property
    def num_parameters(self) -> int:
        """The number of ODE system parameters"""
        return self._pi.num_par

    @property
    def aux_names(self) -> List[str]:
        """The list of auxiliary variable names"""
        return self._pi.aux

    @property
    def num_aux(self) -> int:
        """The number of auxiliary variables"""
        return self._pi.num_aux

    @property
    def num_noise(self) -> int:
        """The number of Wiener variables in the system"""
        return self._pi.num_noise

    @property
    def is_initialized(self) -> bool:
        """Get whether the simulator is initialized"""
        return self._integrator.is_initialized()

    def __init__(
        self,
        variables: Dict[str, float],  # ---ivp---v
        parameters: Dict[str, float],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        src_file: str | None = None,
        rhs_equation: OpenCLRhsEquation | None = None,
        supplementary_equations: List[Callable[[Any], Any]] | None = None,  # ---ivp---^
        single_precision: bool = True,  # --> device/ctx setup?
        stepper: Stepper = Stepper.rk4,  # ---stepper---v
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 1000000,  # ---stepper---^
        max_store: int = 1000000,  # -- trajectory---v
        nout: int = 1,  # -- trajectory---^
        solver_parameters: Optional[SolverParams] = None,
        t_span: Tuple[float, float] = (0.0, 1000.0),
        device_type: CLDeviceType | None = None,  # ---device/ctx setup---
        vendor: CLVendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
        device_ids: List[int] | None = None,
    ) -> None:

        input_file = self._handle_clode_rhs_cl_file(
            src_file, rhs_equation, supplementary_equations
        )

        if aux is None:
            aux = []

        self._pi = ProblemInfo(
            input_file,
            list(variables.keys()),
            list(parameters.keys()),
            aux,
            num_noise,
        )
        self._stepper = stepper
        self._single_precision = single_precision

        # _runtime as an instance variable
        self._runtime = initialize_runtime(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )

        # derived classes override this to call appropriate pybind constructors.
        self._create_integrator()
        self._build_cl_program()

        # set solver_parameters and sync to device
        if solver_parameters is not None:
            self._sp = solver_parameters
        else:
            self._sp = SolverParams(
                dt, dtmax, abstol, reltol, max_steps, max_store, nout
            )
        self.set_solver_parameters()

        # set initial state and parameters, sync to device
        self._variable_defaults = variables
        self._parameter_defaults = parameters
        self._device_initial_state = np.array(
            list(self._variable_defaults.values()), dtype=np.float64, ndmin=2
        )
        self._device_parameters = np.array(
            list(self._parameter_defaults.values()), dtype=np.float64, ndmin=2
        )
        self._ensemble_size = 1
        self._ensemble_shape = (1,)
        self._set_problem_data(self._device_initial_state, self._device_parameters)

        # set tspan and sync to device
        self.set_tspan(t_span=t_span)
        # ---> now the simulator is ready to go

    def _create_integrator(self) -> None:
        self._integrator = SimulatorBase(
            self._pi,
            self._stepper.value,
            self._single_precision,
            self._runtime,
            _clode_root_dir,
        )

    def _build_cl_program(self):
        self._integrator.build_cl()
        self._cl_program_is_valid = True

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

    def set_repeat_ensemble(self, num_repeats: int) -> None:
        """Create an ensemble with identical parameters and initial states.

        This method uses the default parameters and initial state only. For other
        options, see set_ensemble.

        Args:
            num_repeats (int): The number of repeats for the ensemble.

        Returns:
            None
        """

        self._ensemble_size = num_repeats
        self._ensemble_shape = (num_repeats, 1)
        initial_state, parameters = self._make_default_problem_data()
        self._set_problem_data(initial_state=initial_state, parameters=parameters)

    def set_ensemble(
        self,
        variables: Optional[
            Union[np.ndarray, Mapping[str, Union[float, List[float], np.ndarray]]]
        ] = None,
        parameters: Optional[
            Union[np.ndarray, Mapping[str, Union[float, List[float], np.ndarray]]]
        ] = None,
    ) -> None:
        """Set the parameters and/or initial states an ensemble ODE problem, possibly
        changing the ensemble size.

        Generates initial state and parameter arrays with shapes (ensemble_size,
        num_variables) and (ensemble_size, num_parameters), respectively, with one row
        per initial value problem.

        Specifying full arrays or dictionaries mapping parameter/variable names to
        values are supported. The values may be scalars or 1D arrays of a constant
        length. This array length sets the new ensemble_size, and any scalars will be
        broadcast to form fully specified arrays.

        Unspecified values will be taken from the parameter and initial state default
        values. In the case of initial state values, the most recent state from
        simulation will be preferred in the following cases: - when expanding the
        ensemble from size 1 - when the ensemble size does not change

        To override the above behaviour and use the default initial state, specify the
        default initial state as an argument.

        Args:
            variables (np.array | dict): The initial state
            parameters (np.array | dict): The parameters
        """
        self._validate_ensemble_args(variables=variables, parameters=parameters)
        previous_ensemble_size = self._ensemble_size

        var_size = 1
        if isinstance(variables, np.ndarray):
            var_size = variables.shape[0]
        elif isinstance(variables, Mapping):
            vars = {
                self.variables.index(k): v for k, v in parameters.items() if len(v) > 1
            }
            if len(vars) > 0:
                var_size = np.unique([len(v) for v in vars.values()])
                if len(var_size) > 1:
                    raise ValueError(
                        "Arrays specified for initial states must have the same length"
                    )
                var_size = par_size[0]

        par_size = 1
        if isinstance(parameters, np.ndarray):
            par_size = parameters.shape[0]
        elif isinstance(parameters, Mapping):
            pars = {
                self.parameter_names.index(k): v
                for k, v in parameters.items()
                if len(v) > 1
            }
            if len(pars) > 0:
                par_size = np.unique([len(v) for v in pars.values()])
                if len(par_size) > 1:
                    raise ValueError(
                        "Arrays specified for parameters must have the same length"
                    )
                par_size = par_size[0]

        if var_size != 1 and par_size != 1 and var_size != par_size:
            raise ValueError(
                "Arrays specified for parameters and initial states must have the same length"
            )

        self._ensemble_size = max(var_size, par_size)
        self._ensemble_shape = (self._ensemble_size, 1)

        use_current_state = (self._ensemble_size == previous_ensemble_size) | (
            previous_ensemble_size == 1
        )

        vars_array, pars_array = self._make_default_problem_data(
            use_current_state=use_current_state
        )

        if isinstance(parameters, np.ndarray):
            pars_array = parameters
        elif isinstance(parameters, Mapping):
            for index, value in pars.items():
                pars_array[:, index] = np.array(value)

        if isinstance(variables, np.ndarray):
            vars_array = variables
        elif isinstance(variables, Mapping):
            for index, value in vars.items():
                vars_array[:, index] = np.array(value)

        self._set_problem_data(vars_array, pars_array)

    def set_grid_ensemble(
        self,
        variables: Optional[Mapping[str, Union[np.ndarray, float, List[float]]]] = None,
        parameters: Optional[
            Mapping[str, Union[np.ndarray, float, List[float]]]
        ] = None,
        indexing: Literal["xy", "ij"] = "xy",
    ) -> tuple[np.ndarray, ...]:
        """Creates a grid ensmble using Numpy's meshgrid.

        Specify parameters and/or variables with which to make the grid via dictionaries
        mapping names to gridpoint arrays to be used as arguments to meshgrid.

        The resulting ensemble size is the product of lengths of the specified vectors.

        Args:
            variables (mapping[str, array-like], optional): the variables
            parameters (mapping[str, array-like], optional): the parameters
            indexing (['xy', 'ij'], optional): use cartesian ('xy') or matrix ('ij') indexing (see meshgrid)

        Returns:
            None

        """
        # 0. Find ensemble size and shape: shape=(len(p1),len(p2),...), size=prod(shape)
        # 1. generate fully specified arrays with ensemble_size rows
        # 2. replace specified columns with flattened coordinate arrays

        self._validate_ensemble_args(variables=variables, parameters=parameters)

        grid_values = []
        if parameters is not None:
            grid_pars = {
                self.parameter_names.index(k): v
                for k, v in parameters.items()
                if len(v) > 1
            }
            grid_values += list(grid_pars.values())

        if variables is not None:
            grid_vars = {
                self.variable_names.index(k): v
                for k, v in variables.items()
                if len(v) > 1
            }
            grid_values += list(grid_vars.values())

        self._ensemble_shape = tuple([len(v) for v in grid_values])
        self._ensemble_size = np.prod(self._ensemble_shape)

        vars_array, pars_array = self._make_default_problem_data()

        grid_arrays = np.meshgrid(*grid_values, indexing=indexing)

        mesh_index = 0
        if parameters is not None:
            for index in grid_pars.keys():
                pars_array[:, index] = grid_arrays[mesh_index].flatten(order="F")
                mesh_index += 1

        if variables is not None:
            for index in grid_vars.keys():
                vars_array[:, index] = grid_arrays[mesh_index].flatten(order="F")
                mesh_index += 1

        self._set_problem_data(vars_array, pars_array)
        return tuple(grid_arrays)

    def _validate_ensemble_args(
        self,
        variables: Optional[
            Union[np.ndarray, Mapping[str, Union[float, List[float], np.ndarray]]]
        ] = None,
        parameters: Optional[
            Union[np.ndarray, Mapping[str, Union[float, List[float], np.ndarray]]]
        ] = None,
    ) -> None:

        if isinstance(variables, np.ndarray):
            if len(variables.shape) != 2 or variables.shape[1] != self.num_variables:
                raise ValueError(
                    f"initial_state must be a matrix with {self.num_variables} columns"
                )
        elif isinstance(variables, Mapping):
            unknown_variables = set(variables.keys()) - set(self.variable_names)
            if len(unknown_variables) > 0:
                raise ValueError(f"Unknown variable name(s): {unknown_variables}")

        if isinstance(parameters, np.ndarray):
            if len(parameters.shape) != 2 or parameters.shape[1] != self.num_parameters:
                raise ValueError(
                    f"parameters must be a matrix with {self.num_parameters} columns"
                )
        elif isinstance(parameters, Mapping):
            unknown_parameters = set(parameters.keys()) - set(self.parameter_names)
            if len(unknown_parameters) > 0:
                raise ValueError(f"Unknown parameter name(s): {unknown_parameters}")

        if isinstance(variables, np.ndarray) and isinstance(parameters, np.ndarray):
            if variables.shape[0] != parameters.shape[0]:
                raise ValueError(
                    f"initial_state and parameters must have the same number of rows"
                )

    def _make_default_problem_data(
        self, use_current_state: bool = False
    ) -> tuple[np.ndarray, np.ndarray]:
        """Create initial state and parameter arrays from default values

        Args:
            use_current_state (bool, optional): _description_. Defaults to False.

        Returns:
            tuple[np.ndarray, np.ndarray]: the initial state and parameter arrays
        """

        default_initial_state = list(self._variable_defaults.values())
        if use_current_state:
            default_initial_state = self.get_initial_state().squeeze()

        initial_state = np.tile(
            default_initial_state,
            (self._ensemble_size, 1),
        )
        parameters = np.tile(
            list(self._parameter_defaults.values()),
            (self._ensemble_size, 1),
        )
        return initial_state, parameters

    def _set_problem_data(
        self, initial_state: np.ndarray, parameters: np.ndarray
    ) -> None:
        """Set both initial state and parameters at the same time.

        This method supports changing ensemble size, but initial state and parameters
        must be completely specified as ndarrays with the same number of rows.

        Args:
            initial_state (np.array): The initial state. shape=(ensemble_size, num_variables)
            parameters (np.array): The parameters. shape=(ensemble_size, num_parameters)
        """

        self._device_initial_state = initial_state
        self._device_parameters = parameters
        self._integrator.set_problem_data(
            initial_state.flatten(order="F"),
            parameters.flatten(order="F"),
        )

    def _set_parameters(self, parameters: np.ndarray) -> None:
        """Set the ensemble parameters without changing ensemble size.

        New ensemble parameters must match the current ensemble size.

        Args:
            parameters (np.array): The parameters. shape=(ensemble_size, num_parameters)
        """
        self._device_parameters = parameters
        self._integrator.set_pars(parameters.flatten(order="F"))

    def _set_initial_state(self, initial_state: np.ndarray) -> None:
        """Set the initial state without changing ensemble size.

        New ensemble initial_state must match the current ensemble size.

        Args:
            initial_state (np.ndarray): The initial state. shape=(ensemble_size, num_variables)
        """
        self._device_initial_state = initial_state
        self._integrator.set_x0(initial_state.flatten(order="F"))

    def set_tspan(self, t_span: tuple[float, float]) -> None:
        """Set the time span of the simulation.

        Args:
            t_span (tuple[float, float]): The time span.
        """
        self._t_span = t_span
        self._integrator.set_tspan(t_span)

    def get_tspan(self) -> tuple[float, float]:
        """Returns the simulation time span currently set on the device.

        Returns:
            tuple[float, float]: The time span
        """
        self._t_span = self._integrator.get_tspan()

    def shift_tspan(self) -> None:
        """Shift the time span to the current time plus the time period."""
        self._integrator.shift_tspan()
        self._t_span = self._integrator.get_tspan()

    def set_solver_parameters(
        self,
        solver_parameters: Optional[SolverParams] = None,
        dt: Optional[float] = None,
        dtmax: Optional[float] = None,
        abstol: Optional[float] = None,
        reltol: Optional[float] = None,
        max_steps: Optional[int] = None,
        max_store: Optional[int] = None,
        nout: Optional[int] = None,
    ) -> None:
        """Update solver parameters and push to the device.

        A full solver parameters struct or individual fields may be specified

        Args:
            solver_parameters (SolverParams, optional): A solver parameters structure. Defaults to None.
            dt (float, optional): The time step. Defaults to None.
            dtmax (float, optional): Maximum time step for adaptive solvers. Defaults to None.
            abstol (float, optional): Absolute tolerance for adaptive solvers. Defaults to None.
            reltol (float, optional): Relative tolerance for adaptive solvers. Defaults to None.
            max_steps (int, optional): Maximum number of time steps. Defaults to None.
            max_store (int, optional): Maximum steps to store for trajectories. Defaults to None.
            nout (int, optional): Store interval, in number of steps, for trajectories. Defaults to None.
        """
        if solver_parameters is not None:
            self._sp = solver_parameters
        else:
            if dt is not None:
                self._sp.dt = dt
            if dtmax is not None:
                self._sp.dtmax = dtmax
            if abstol is not None:
                self._sp.abstol = abstol
            if reltol is not None:
                self._sp.reltol = reltol
            if max_steps is not None:
                self._sp.max_steps = max_steps
            if max_store is not None:
                self._sp.max_store = max_store
            if nout is not None:
                self._sp.nout = nout
        self._integrator.set_solver_params(self._sp)

    def get_solver_parameters(self):
        """Get the current ensemble parameters from the OpenCL device

        Returns:
            SolverParams: The solver parameters structure
        """
        return self._integrator.get_solver_params()

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

    def transient(
        self,
        t_span: Optional[Tuple[float, float]] = None,
        update_x0: bool = True,
        fetch_results: bool = False,
    ) -> Optional[np.ndarray]:
        """Run a transient simulation.

        Args:
            t_span (tuple[float, float]): Time interval for integration.
            update_x0 (bool, optional): Whether to update the initial state. Defaults to True.
            fetch_results (bool): Whether to fetch the feature results from the device and return them here

        Returns:
            None
        """

        # Lazy rebuild - would also need to verify device data is set
        # if not self._cl_program_is_valid:
        #     self._integrator.build_cl()
        #     self._cl_program_is_valid = True

        # this should be a verification that all args have been set (lazily do so here if needed?)
        # if not self.is_initialized:
        #     raise RuntimeError("Simulator is not initialized")

        if t_span is not None:
            self.set_tspan(t_span=t_span)

        self._integrator.transient()
        # invalidates _device_final_state
        self._device_final_state = None

        if update_x0:
            self._integrator.shift_x0()
            # invalidates _device_initial_state (?)
            self._device_initial_state = None

        if fetch_results:
            return self.get_final_state()
            # Note that this triggers a device-to-host transfer.

    def get_initial_state(self) -> np.ndarray:
        """Get the initial state of the simulation from the device.

        Note that this triggers a device-to-host transfer.

        Returns:
            np.array: The initial state of the simulation.
        """
        if self._device_initial_state is None:
            initial_state = np.array(self._integrator.get_x0(), dtype=np.float64)
            self._device_initial_state = initial_state.reshape(
                (self._ensemble_size, self.num_variables), order="F"
            )
        return self._device_initial_state

    def get_final_state(self) -> np.ndarray:
        """Get the final state of the simulation from the device.

        Note that this triggers a device-to-host transfer.

        Returns:
            np.array: The final state of the simulation.
        """
        if self._device_final_state is None:
            final_state = np.array(self._integrator.get_xf(), dtype=np.float64)
            self._device_final_state = final_state.reshape(
                (self._ensemble_size, self.num_variables), order="F"
            )
        return self._device_final_state

    def get_max_memory_alloc_size(self, deviceID: int = 0) -> int:
        """Get the device maximum memory allocation size

        Args:
            deviceID (int, optional): The device ID. Defaults to 0.

        Returns:
            int: The maximum size of memory allocation in GB
        """
        return self._runtime.get_max_memory_alloc_size(deviceID)

    def get_double_support(self, deviceID: int = 0) -> bool:
        """Get whether the device supports double precision

        Args:
            deviceID (int, optional): The device ID. Defaults to 0.

        Returns:
            bool: Whether the device supports double precision
        """
        return self._runtime.get_double_support(deviceID)

    def get_device_cl_version(self, deviceID: int = 0) -> str:
        """Get the device OpenCL version

        Args:
            deviceID (int, optional): The device ID. Defaults to 0.

        Returns:
            str: the device CL version
        """
        return self._runtime.get_device_cl_version(deviceID)

    def get_available_steppers(self) -> List[str]:
        """Get the list of valid time stepper names"""
        return self._integrator.get_available_steppers()

    def get_program_string(self) -> str:
        """Get the clODE OpenCL program string"""
        return self._integrator.get_program_string()

    def print_status(self) -> None:
        """Print the simulator status info"""
        old_level = get_log_level()
        set_log_level(LogLevel.info)
        self._integrator.print_status()
        set_log_level(old_level)

    def print_devices(self) -> None:
        """Print the available devices"""
        old_level = get_log_level()
        set_log_level(LogLevel.info)
        self._runtime.print_devices()
        set_log_level(old_level)
