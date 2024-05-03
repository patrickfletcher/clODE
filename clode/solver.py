from __future__ import annotations

from enum import Enum
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np

from .function_converter import OpenCLConverter, OpenCLRhsEquation
from .runtime import (
    CLDeviceType,
    CLVendor,
    OpenCLResource,
    _clode_root_dir,
    initialize_runtime,
)
from clode.cpp.clode_cpp_wrapper import SolverParams, ProblemInfo, SimulatorBase

from .xpp_parser import convert_xpp_file


# TODO[API]: different steppers use different parameters subsets. Expose only relevant. Model with classes/mixins?
# - StepperBase + mixins --> each stepper as python class. Dict {"name": class} from "get_available_steppers()"
class Stepper(Enum):
    euler = "euler"
    heun = "heun"
    rk4 = "rk4"
    bs23 = "bs23"
    dormand_prince = "dopri5"
    stochastic_euler = "seuler"

# TODO[API]: defaults for solverParams are copied in each constructor plus wrapper. Should be only ONE place globally.
# - Prefer the struct defaults? Create a default config?
DEFAULT_SOLVER_PARAMS = SolverParams(dt=0.1, dtmax=1.0, abstol=1e-6, reltol=1e-4, max_steps=100000, max_store=100000, nout=1)

class Simulator:
    """Base class for simulating an ensemble of instances of an ODE system.

    It provides the core functionality for advancing the simulation in time without storing any intermediate state.  
    It can be used directly if only the final state is of interest, or as a base class for other simulators.
    """

    _integrator: SimulatorBase
    _runtime: OpenCLResource
    _single_precision: bool
    _stepper: Stepper
    _pi: ProblemInfo

    # changes to the above items require rebuilding the CL program. Flag them and rebuild if necessary on simulation function call
    _cl_program_is_valid: bool = False

    # changes to the remaining items do not require rebuilding. Setters ensure device memory is correctly allocated and populated with valid data
    _sp: SolverParams
    _t_span: Tuple[float, float]

    # values may change, but not keys:
    _variable_defaults: Dict[str, float]
    _parameter_defaults: Dict[str, float]

    _variables: Optional[Dict[str, np.ndarray]] = None 
    _parameters: Optional[Dict[str, np.ndarray]] = None

    _ensemble_size: int # C++ layer: nPts
    _ensemble_shape: Tuple #support 1D (anything), and 2D+ grids. Outputs can be reshaped to match (e.g., 2D grid on two parameters, return state, features as 2D)
    _ensemble_parameters: Optional[np.ndarray] = None     #2D array shape (ensemble_size, num_parameters)
    _ensemble_initial_state: Optional[np.ndarray] = None  #2D array shape (ensemble_size, num_variables)
    _ensemble_final_state: Optional[np.ndarray] = None    #2D array shape (ensemble_size, num_variables)

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
        variables: Dict[str, float], #---ivp---
        parameters: Dict[str, float],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,         
        src_file: str | None = None,
        rhs_equation: OpenCLRhsEquation | None = None,
        supplementary_equations: List[Callable[[Any], Any]] | None = None, #---ivp---    
        single_precision: bool = True,  #-> with device???
        stepper: Stepper = Stepper.rk4, #---stepper---encapsulate with external class, expose only relevant params
        dt: float = 0.1, 
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 1000000,       #---stepper---
        max_store: int = 1000000,   # -- only for trajectory
        nout: int = 1,              # -- only for trajectory
        solver_parameters: SolverParams = None,
        t_span: Tuple[float, float] = (0.0, 1000.0),  # more natural to set in simulation function calls
        device_type: CLDeviceType | None = None, #---device/ctx setup--- feels annoying here [do like pyopencl: "make_some_context()"]
        vendor: CLVendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
        device_ids: List[int] | None = None,
    ) -> None:
        
        input_file = self._handle_clode_rhs_cl_file(
            src_file, rhs_equation, supplementary_equations
        )

        self._stepper = stepper
        self._single_precision = single_precision
        self._sp = SolverParams(dt, dtmax, abstol, reltol, max_steps, max_store, nout)
        self._t_span = t_span

        self._variable_defaults = variables
        self._parameter_defaults = parameters
        self._ensemble_size = 1
        self._ensemble_shape = (1,)

        if aux is None:
            aux = []

        self._pi = ProblemInfo(
            input_file,
            list(variables.keys()),
            list(parameters.keys()),
            aux,
            num_noise,
        )

        # _runtime as an instance variable
        self._runtime = initialize_runtime(
            device_type,
            vendor,
            platform_id,
            device_id,
            device_ids,
        )

        self._build_integrator()
        self._init_integrator()


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
    

    # Ensemble setter methods

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
        self._ensemble_shape = (num_repeats,)
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

        self._ensemble_size = cl_array_length
        self._variables = self._create_cl_arrays(local_variables, cl_array_length)
        self._parameters = self._create_cl_arrays(local_parameters, cl_array_length)

        # if self.is_initialized:
        #     vars_array, pars_array = self._pack_data()
        #     self._integrator.set_problem_data(vars_array, pars_array)
        # else:
        self._init_integrator()
        # self.seed_rng(seed)


    # TODO: support passing arrays with shapes, and return matching shape for results?
    # - e.g., two 2D arrays from meshgrid for two-parameter sweep
    def set_grid_ensemble(self, parameters, variables):
        ''' Creates a grid ensmble using Numpy's meshgrid.

        '''
        # np.meshgrid
        pass


    def _validate_problem_data(self, parameters, variables):
        pass


    def _build_integrator(self) -> None:
        self._integrator = SimulatorBase(
            self._pi,
            self._stepper.value,
            self._single_precision,
            self._runtime,
            _clode_root_dir,
        )
        self._integrator.build_cl()
        self._cl_program_is_valid = True

    def _init_integrator(self) -> None:
        # FeaturesSimulator overrides this
        vars_array, pars_array = self._pack_data()
        self._integrator.initialize(
            self._t_span,
            vars_array,
            pars_array,
            self._sp,
        )


    # Low-level functions for granular control.
    # - verify that new value array matches (or can be broadcast to match) current ensemble_size
    # - ensure correct order/packing (as in pack_data), flatten. 

    def set_problem_data(self, initial_state: np.array, parameters: np.array) -> None:
        """Set the problem data.

        Set both initial state and parameters at the same time. 
        This method supports changing ensemble size and does not require re-building CL program

        Args:
            initial_state (np.array): The initial state.
            parameters (np.array): The parameters.

        Returns:
            None
        """
        
        raise NotImplementedError

    def set_parameters(self, parameters: np.ndarray) -> None:
        """Set the parameters.

        New ensemble parameters must match the current initial state dimensions.

        Parameters
        ----------
            parameters (np.array): The parameters.

        Returns:
        ----------
            None
        """
        # self._device_parameters = parameters
        # self._integrator.set_pars(parameters)
        pass


    def set_initial_state(
        self, initial_state: dict[str, Union[float, np.ndarray, List[float]]]
    ) -> None:
        """Set the initial state.

        New ensemble initial_state must match the current ensemble parameter dimensions.

        Args:
            initial_state (dict[str, Union[float, np.ndarray, List[float]]]): The initial state.

        Returns:
            None
        """
        # self._device_initial_state = initial_state
        pass


    # Changing tspan does not require re-building CL program
    def set_tspan(self, t_span: Tuple[float, float]) -> None:
        """Set the time span to simulate over.

        Args:
            t_span (Tuple[float, float]): The time span to simulate over.

        Returns:
            None
        """
        self._t_span = t_span
        self._integrator.set_tspan(t_span)

    def get_tspan(self) -> tuple[float,float]:
        '''Returns the simulation time span currently set on the device'''
        self._t_span = self._integrator.get_tspan()
    

    def shift_tspan(self) -> None:
        """Shift the time span to the current time plus the time period."""
        self._integrator.shift_tspan()
        self._t_span = self._integrator.get_tspan()


    # Changing solver parameters does not require re-building CL program
    # NOTE: int inputs to float are autocast, but not float --> int. Should we cast here and warn?
    # TODO: would be nice to separate options based on stepper types. 
    # - abstol/reltol/dtmax only relevant to adaptive solvers. Additional options will be needed when steppers are expanded...
    # - Similarly, max_store/nout are only relevant to trajectory
    def set_solver_parameters(
        self, 
        sp: Optional[SolverParams] = None,
        dt: Optional[float] = None, 
        dtmax: Optional[float] = None, 
        abstol: Optional[float] = None, 
        reltol: Optional[float] = None, 
        max_steps: Optional[int] = None, 
        max_store: Optional[int] = None, 
        nout: Optional[int] = None,
    ) -> None:
        """Update any of the solver parameters and push to device

        Args:

        Returns:
            None
        """
        if sp is not None:
            self._sp = sp
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
        return self._integrator.get_solver_params()
    
    # Seeding RNG does not require re-building CL program
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
    

    def transient(self, 
                  t_span:Optional[Tuple[float, float]] = None, 
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

        #this should be a verification that all args have been set (lazily do so here if needed?)
        if not self.is_initialized:
            raise RuntimeError("Simulator is not initialized")
        
        if t_span is not None:
            self.set_tspan(t_span=t_span)

        self._integrator.transient()
        # invalidates _ensemble_final_state
        self._ensemble_final_state = None

        if update_x0:
            self._integrator.shift_x0()
            # invalidates _ensemble_initial_state (?)
            self._ensemble_initial_state = None

        if fetch_results:
            return self.get_final_state()
            # Note that this triggers a device-to-host transfer.



    def get_initial_state(self) -> np.ndarray:
        """Get the initial state of the simulation from the device.

        Note that this triggers a device-to-host transfer.

        Returns:
            np.array: The initial state of the simulation.
        """
        if self._ensemble_initial_state is None:    
            initial_state = np.array(self._integrator.get_x0(), dtype=np.float64)
            self._ensemble_initial_state = initial_state.reshape(
                (self._ensemble_size, self.num_variables),
                order='F')
        return self._ensemble_initial_state


    def get_final_state(self) -> np.ndarray:
        """Get the final state of the simulation from the device.

        Note that this triggers a device-to-host transfer.

        Returns:
            np.array: The final state of the simulation.
        """
        if self._ensemble_final_state is None:    
            final_state = np.array(self._integrator.get_xf(), dtype=np.float64)
            self._ensemble_final_state = final_state.reshape(
                (self._ensemble_size, self.num_variables),
                order='F')
        return self._ensemble_final_state



    def get_max_memory_alloc_size(self) -> int:
        """Get the device maximum memory allocation size"""
        return self._runtime.get_max_memory_alloc_size()

    def get_double_support(self) -> bool:
        """Get whether the device supports double precision"""
        return self._runtime.get_double_support()
    

    def get_available_steppers(self) -> List[str]:
        """Get the list of valid time stepper names"""
        return self._integrator.get_available_steppers()

    def get_device_cl_version(self) -> str:
        """Get the device OpenCL version"""
        return self._runtime.get_device_cl_version()
    
    def get_program_string(self) -> str:
        """Get the clODE OpenCL program string"""
        return self._integrator.get_program_string()

    #TODO: this requires correct logger setting - default doesn't print
    def print_status(self) -> None:
        """Print the simulator status info"""
        self._integrator.print_status()

    #TODO: this requires correct logger setting - default doesn't print
    def print_devices(self) -> None:
        """Print the available devices"""
        self._runtime.print_devices()