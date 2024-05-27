from __future__ import annotations

from enum import Enum

# TODO[typing] - standardize typing throughout the package
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np

# - npt.NDArray[np.float64], npt.ArrayLike
# https://numpy.org/neps/nep-0029-deprecation_policy.html
import numpy.typing as npt

from clode.cpp.clode_cpp_wrapper import ProblemInfo, SimulatorBase, SolverParams

from .function_converter import OpenCLConverter, OpenCLRhsEquation
from .runtime import (
    CLDeviceType,
    CLVendor,
    LogLevel,
    OpenCLResource,
    _clode_root_dir,
    get_log_level,
    initialize_runtime,
    set_log_level,
)
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
# TODO[API]: set_ensemble and set_repeats
# when possible, convention to use (keep/broadcast) previous params/initial_state?
# - goes with convention to advance initial_state by default
# - when not possible, the defaults will be used as basis for new ensembles
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

    _ensemble_size: int  # C++ layer: nPts
    _ensemble_shape: Tuple

    # 2D array shape (ensemble_size, num_parameters)
    _device_parameters: Optional[np.ndarray] = None

    # 2D array shape (ensemble_size, num_variables)
    _device_initial_state: Optional[np.ndarray] = None
    _device_final_state: Optional[np.ndarray] = None

    # TODO[mkdocs] - put these below init?
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

    def __init__(
        self,
        variables: Dict[str, float],
        parameters: Dict[str, float],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        src_file: Optional[str] = None,
        rhs_equation: Optional[OpenCLRhsEquation] = None,
        supplementary_equations: List[Callable[[Any], Any]] | None = None,
        stepper: Stepper = Stepper.rk4,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 1000000,
        max_store: int = 1000000,
        nout: int = 1,
        solver_parameters: Optional[SolverParams] = None,
        t_span: Tuple[float, float] = (0.0, 1000.0),
        single_precision: bool = True,
        device_type: Optional[CLDeviceType] = None,
        vendor: Optional[CLVendor] = None,
        platform_id: Optional[int] = None,
        device_id: Optional[int] = None,
        device_ids: Optional[List[int]] = None,
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

        # set tspan and sync to device
        self.set_tspan(t_span=t_span)

        # set initial state and parameters, sync to device
        self._variable_defaults = variables
        self._parameter_defaults = parameters

        # use set_repeat_ensemble(1)?
        self._ensemble_size = 1
        self._ensemble_shape = (1,)
        default_initial_state = np.array(
            list(self._variable_defaults.values()), dtype=np.float64, ndmin=2
        )
        default_parameters = np.array(
            list(self._parameter_defaults.values()), dtype=np.float64, ndmin=2
        )
        self._set_problem_data(default_initial_state, default_parameters)

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
        initial_state, parameters = self._make_problem_data(
            new_size=num_repeats, new_shape=(num_repeats, 1)
        )
        self._set_problem_data(initial_state=initial_state, parameters=parameters)

    # TODO: refactor some parts?
    # TODO[typing]
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
        if variables is None and parameters is None:
            raise ValueError(f"initial_state and parameters cannot both be None")

        # validate variables argument
        if isinstance(variables, np.ndarray):
            if len(variables.shape) != 2 or variables.shape[1] != self.num_variables:
                raise ValueError(
                    f"initial_state must be a matrix with {self.num_variables} columns"
                )
        elif isinstance(variables, Mapping):
            unknown_variables = set(variables.keys()) - set(self.variable_names)
            if len(unknown_variables) > 0:
                raise ValueError(f"Unknown variable name(s): {unknown_variables}")
        elif variables is not None:
            raise ValueError(
                f"Expected np.ndarray or Mapping for variables, but got {type(variables)}"
            )

        # validate parameters argument
        if isinstance(parameters, np.ndarray):
            if len(parameters.shape) != 2 or parameters.shape[1] != self.num_parameters:
                raise ValueError(
                    f"parameters must be a matrix with {self.num_parameters} columns"
                )
        elif isinstance(parameters, Mapping):
            unknown_parameters = set(parameters.keys()) - set(self.parameter_names)
            if len(unknown_parameters) > 0:
                raise ValueError(f"Unknown parameter name(s): {unknown_parameters}")
        elif parameters is not None:
            raise ValueError(
                f"Expected np.ndarray or Mapping for parameters, but got {type(variables)}"
            )

        # get the shape and size from variables
        var_size = 1
        var_shape = (1,)
        if isinstance(variables, np.ndarray):
            var_size = variables.shape[0]
            var_shape = (var_size, 1)
        elif isinstance(variables, Mapping):
            variables = {k: np.array(v, dtype=np.float64) for k, v in variables.items()}
            # size/shape from dict. scalars have size=1, shape=()
            var_sizes = [v.size for v in variables.values() if v.size > 1]
            var_shapes = [v.shape for v in variables.values() if v.size > 1]
            if len(set(var_shapes)) > 1:
                shapes = {k: v.shape for k, v in variables.items() if v.size > 1}
                raise ValueError(f"Shape of arrays for variables don't match: {shapes}")
            if len(var_sizes) > 0:
                var_size = var_sizes[0]
                var_shape = var_shapes[0]

        # get the shape and size from parameters
        par_size = 1
        par_shape = (1,)
        if isinstance(parameters, np.ndarray):
            par_size = parameters.shape[0]
            par_shape = (par_size, 1)
        elif isinstance(parameters, Mapping):
            parameters = {
                k: np.array(v, dtype=np.float64) for k, v in parameters.items()
            }
            # size/shape from dict. scalars have size=1, shape=()
            par_sizes = [v.size for v in parameters.values() if v.size > 1]
            par_shapes = [v.shape for v in parameters.values() if v.size > 1]
            if len(set(par_shapes)) > 1:
                shapes = {k: v.shape for k, v in parameters.items() if v.size > 1}
                raise ValueError(
                    f"Shape of arrays for parameters don't match: {shapes}"
                )
            if len(par_sizes) > 0:
                par_size = par_sizes[0]
                par_shape = par_shapes[0]

        # size and shape must match
        # TODO: shape check? only for dict case...
        if var_size > 1 and par_size > 1:
            if var_size != par_size or var_size != par_size:
                raise ValueError(
                    "Arrays specified for parameters and initial states must have the same size"
                )

        # print(var_size, var_shape, par_size, par_shape)
        new_size = var_size if var_size > 1 else par_size
        new_shape = var_shape if var_size > 1 else par_shape

        vars_array, pars_array = self._make_problem_data(
            variables=variables,
            parameters=parameters,
            new_size=new_size,
            new_shape=new_shape,
        )
        self._set_problem_data(vars_array, pars_array)

    # TODO: when to keep/broadcast current vs default values?
    # TODO[typing]
    def _make_problem_data(
        self,
        variables: Optional[dict[str, np.ndarray]] = None,
        parameters: Optional[dict[str, np.ndarray]] = None,
        new_size: Optional[int] = None,
        new_shape: Optional[tuple[int, ...]] = None,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Create initial state and parameter arrays from default values

        The resulting arrays by convention have shapes (ensemble_size, num_variables)
        and (ensemble_size, num_parameters)

        Args:
            use_current_state (bool, optional): _description_. Defaults to False.

        Returns:
            tuple[np.ndarray, np.ndarray]: the initial state and parameter arrays
        """

        if len(new_shape) == 1:
            new_shape = (new_size, 1)

        previous_size = self._ensemble_size
        valid_previous_size = (previous_size == new_size) | (previous_size == 1)

        if valid_previous_size:
            initial_state_array = self.get_initial_state()
            parameter_array = self._device_parameters
        else:
            initial_state_array = np.array(
                list(self._variable_defaults.values()), dtype=np.float64, ndmin=2
            )
            parameter_array = np.array(
                list(self._parameter_defaults.values()), dtype=np.float64, ndmin=2
            )

        if initial_state_array.shape[0] == 1:
            initial_state_array = np.tile(initial_state_array, new_shape)

        if parameter_array.shape[0] == 1:
            parameter_array = np.tile(parameter_array, new_shape)

        # possibly overwrite some or all of the arrays
        if isinstance(variables, np.ndarray):
            initial_state_array = variables
        elif isinstance(variables, Mapping):
            for key, value in variables.items():
                index = self.variable_names.index(key)
                value = np.repeat(value, new_size) if value.size == 1 else value
                initial_state_array[:, index] = np.array(value.flatten())

        if isinstance(parameters, np.ndarray):
            parameter_array = parameters
        elif isinstance(parameters, Mapping):
            for key, value in parameters.items():
                index = self.parameter_names.index(key)
                value = np.repeat(value, new_size) if value.size == 1 else value
                parameter_array[:, index] = np.array(value.flatten())

        self._ensemble_size = new_size
        self._ensemble_shape = new_shape
        return initial_state_array, parameter_array

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
        # old_level = get_log_level()
        # set_log_level(LogLevel.info)
        self._integrator.print_status()
        # set_log_level(old_level)

    def print_devices(self) -> None:
        """Print the available devices"""
        # old_level = get_log_level()
        # set_log_level(LogLevel.info)
        self._runtime.print_devices()
        # set_log_level(old_level)
