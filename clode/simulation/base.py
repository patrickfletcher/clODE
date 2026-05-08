from __future__ import annotations

from enum import Enum
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt

from .._opencl.factory import create_simulator_backend
from ..problem.definition import ProblemInfo
from ..problem.python import OpenCLConverter, OpenCLRhsEquation
from ..problem.source import RhsSource, create_rhs_source, load_rhs_source
from ..problem.xpp import convert_xpp_file
from ..runtime import (
	CLDeviceType,
	CLVendor,
	LogLevel,
	OpenCLResource,
	_clode_root_dir,
	get_log_level,
	initialize_runtime,
	set_log_level,
)
from ..runtime.selection import RuntimeSelection
from ._protocols import SimulatorBackend
from .params import SolverParams


class Stepper(Enum):
	euler = "euler"
	heun = "heun"
	rk4 = "rk4"
	bs23 = "bs23"
	dormand_prince = "dopri5"
	stochastic_euler = "seuler"


class Simulator:
	"""Base class for simulating an ensemble of instances of an ODE system.

	It provides the core functionality for advancing the simulation in time without
	storing any intermediate state. May be used directly when only the final state is
	of interest, or as a base class for other simulators.
	"""

	_integrator: SimulatorBackend
	_runtime: OpenCLResource
	_runtime_selection: RuntimeSelection
	_single_precision: bool
	_stepper: Stepper
	_pi: ProblemInfo
	_rhs_source: RhsSource

	_cl_program_is_valid: bool = False

	_sp: SolverParams
	_t_span: Tuple[float, float]

	_variable_defaults: Dict[str, float]
	_parameter_defaults: Dict[str, float]

	_ensemble_size: int
	_ensemble_shape: Tuple

	_device_parameters: Optional[np.ndarray] = None
	_device_initial_state: Optional[np.ndarray] = None
	_device_final_state: Optional[np.ndarray] = None
	_device_dt: Optional[np.ndarray] = None
	_device_tf: Optional[np.ndarray] = None

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

		self._rhs_source = self._prepare_rhs_source(
			src_file, rhs_equation, supplementary_equations
		)

		if aux is None:
			aux = []

		self._pi = ProblemInfo(
			self._rhs_source.origin_label,
			list(variables.keys()),
			list(parameters.keys()),
			aux,
			num_noise,
		)
		self._stepper = stepper
		self._single_precision = single_precision
		self._runtime_selection = RuntimeSelection(
			device_type=device_type,
			vendor=vendor,
			platform_id=platform_id,
			device_id=device_id,
			device_ids=None if device_ids is None else tuple(device_ids),
		)

		self._runtime = initialize_runtime(
			device_type,
			vendor,
			platform_id,
			device_id,
			device_ids,
		)

		self._create_integrator()
		self._build_cl_program()

		if solver_parameters is not None:
			self._sp = solver_parameters
		else:
			self._sp = SolverParams(
				dt, dtmax, abstol, reltol, max_steps, max_store, nout
			)
		self.set_solver_parameters()

		self.set_tspan(t_span=t_span)

		self._variable_defaults = variables
		self._parameter_defaults = parameters

		self._ensemble_size = 1
		self._ensemble_shape = (1,)
		default_initial_state = np.array(
			list(self._variable_defaults.values()), dtype=np.float64, ndmin=2
		)
		default_parameters = np.array(
			list(self._parameter_defaults.values()), dtype=np.float64, ndmin=2
		)
		self._set_problem_data(default_initial_state, default_parameters)

	def _create_integrator(self) -> None:
		self._integrator = create_simulator_backend(
			self._pi,
			self._rhs_source,
			self._stepper.value,
			self._single_precision,
			_clode_root_dir,
			runtime_selection=self._runtime_selection,
		)

	def _build_cl_program(self):
		self._integrator.build_cl()
		self._cl_program_is_valid = True

	def _ensure_cl_program(self) -> None:
		if not self._cl_program_is_valid:
			self._build_cl_program()

	def _invalidate_solution_cache(self) -> None:
		self._device_final_state = self._device_dt = self._device_tf = None

	def _prepare_rhs_source(
		self,
		src_file: str | None = None,
		rhs_equation: OpenCLRhsEquation | None = None,
		supplementary_equations: List[Callable[[Any], Any]] | None = None,
	) -> RhsSource:

		if src_file is not None and rhs_equation is not None:
			raise ValueError("Cannot specify both src_file and rhs_equation")
		elif src_file is not None:
			if src_file.endswith(".xpp"):
				input_file = convert_xpp_file(src_file)
			else:
				input_file = src_file
			return load_rhs_source(input_file)
		elif rhs_equation is not None:
			converter = OpenCLConverter()
			if supplementary_equations is not None:
				for eq in supplementary_equations:
					converter.convert_to_opencl(eq)
			eqn = converter.convert_to_opencl(
				rhs_equation, mutable_args=[3, 4], function_name="getRHS"
			)
			return create_rhs_source("clode_rhs.cl", eqn)
		else:
			raise ValueError("Must specify either src_file or rhs_equation")

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

		var_size = 1
		var_shape = (1,)
		if isinstance(variables, np.ndarray):
			var_size = variables.shape[0]
			var_shape = (var_size, 1)
		elif isinstance(variables, Mapping):
			variables = {k: np.array(v, dtype=np.float64) for k, v in variables.items()}
			var_sizes = [v.size for v in variables.values() if v.size > 1]
			var_shapes = [v.shape for v in variables.values() if v.size > 1]
			if len(set(var_shapes)) > 1:
				shapes = {k: v.shape for k, v in variables.items() if v.size > 1}
				raise ValueError(f"Shape of arrays for variables don't match: {shapes}")
			if len(var_sizes) > 0:
				var_size = var_sizes[0]
				var_shape = var_shapes[0]

		par_size = 1
		par_shape = (1,)
		if isinstance(parameters, np.ndarray):
			par_size = parameters.shape[0]
			par_shape = (par_size, 1)
		elif isinstance(parameters, Mapping):
			parameters = {
				k: np.array(v, dtype=np.float64) for k, v in parameters.items()
			}
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

		if var_size > 1 and par_size > 1:
			if var_size != par_size or var_size != par_size:
				raise ValueError(
					"Arrays specified for parameters and initial states must have the same size"
				)

		new_size = var_size if var_size > 1 else par_size
		new_shape = var_shape if var_size > 1 else par_shape

		vars_array, pars_array = self._make_problem_data(
			variables=variables,
			parameters=parameters,
			new_size=new_size,
			new_shape=new_shape,
		)
		self._set_problem_data(vars_array, pars_array)

	def _make_problem_data(
		self,
		variables: Optional[dict[str, np.ndarray]] = None,
		parameters: Optional[dict[str, np.ndarray]] = None,
		new_size: Optional[int] = None,
		new_shape: Optional[tuple[int, ...]] = None,
	) -> tuple[np.ndarray, np.ndarray]:
		"""Create initial state and parameter arrays from default values.

		The resulting arrays by convention have shapes (ensemble_size, num_variables)
		and (ensemble_size, num_parameters)
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
			initial_state_array = np.tile(initial_state_array, (new_size, 1))

		if parameter_array.shape[0] == 1:
			parameter_array = np.tile(parameter_array, (new_size, 1))

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
		"""Set both initial state and parameters at the same time."""
		self._device_initial_state = initial_state
		self._device_parameters = parameters
		self._integrator.set_problem_data(
			initial_state.flatten(order="F"),
			parameters.flatten(order="F"),
		)

	def _set_parameters(self, parameters: np.ndarray) -> None:
		"""Set the ensemble parameters without changing ensemble size."""
		self._device_parameters = parameters
		self._integrator.set_pars(parameters.flatten(order="F"))

	def _set_initial_state(self, initial_state: np.ndarray) -> None:
		"""Set the initial state without changing ensemble size."""
		self._device_initial_state = initial_state
		self._integrator.set_x0(initial_state.flatten(order="F"))

	def set_tspan(self, t_span: tuple[float, float]) -> None:
		"""Set the time span of the simulation."""
		self._t_span = t_span
		self._integrator.set_tspan(t_span)

	def get_tspan(self) -> tuple[float, float]:
		"""Returns the simulation time span currently set on the device."""
		self._t_span = tuple(self._integrator.get_tspan())
		return self._t_span

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
		"""Update solver parameters and push to the device."""
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
		self._device_dt = None

	def get_solver_parameters(self) -> SolverParams:
		"""Get the current ensemble parameters from the OpenCL device."""
		return self._integrator.get_solver_params()

	def seed_rng(self, seed: int | None = None) -> None:
		"""Seed the random number generator."""

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
		"""Run a transient simulation."""

		self._ensure_cl_program()

		if t_span is not None:
			self.set_tspan(t_span=t_span)

		self._integrator.transient()
		self._invalidate_solution_cache()

		if update_x0:
			self._integrator.shift_x0()
			self._device_initial_state = None

		if fetch_results:
			return self.get_final_state()

	def get_initial_state(self) -> np.ndarray:
		"""Get the initial state of the simulation from the device."""
		if self._device_initial_state is None:
			self._device_initial_state = np.array(
				self._integrator.get_x0(), dtype=np.float64
			).reshape((self._ensemble_size, self.num_variables), order="F")
		return self._device_initial_state

	def get_final_state(self) -> np.ndarray:
		"""Get the final state of the simulation from the device."""
		if self._device_final_state is None:
			final_state = self._integrator.get_xf()

			if final_state is None:
				raise ValueError("Must run a simulation before getting final state")

			self._device_final_state = np.array(final_state, dtype=np.float64).reshape(
				(self._ensemble_size, self.num_variables), order="F"
			)
		return self._device_final_state

	def get_dt(self) -> np.ndarray:
		"""Get the array of timestep sizes (dt) from the device."""
		if self._device_dt is None:
			self._device_dt = np.array(
				self._integrator.get_dt(), dtype=np.float64
			).reshape(self._ensemble_shape, order="F")
		return self._device_dt

	def get_final_time(self) -> np.ndarray:
		"""Get the array of final times from the device."""
		if self._device_tf is None:
			self._device_tf = np.array(
				self._integrator.get_tf(), dtype=np.float64
			).reshape(self._ensemble_shape, order="F")
		return self._device_tf

	def get_max_memory_alloc_size(self, deviceID: int = 0) -> int:
		"""Get the device maximum memory allocation size."""
		return self._runtime.get_max_memory_alloc_size(deviceID)

	def get_double_support(self, deviceID: int = 0) -> bool:
		"""Get whether the device supports double precision."""
		return self._runtime.get_double_support(deviceID)

	def get_device_cl_version(self, deviceID: int = 0) -> str:
		"""Get the device OpenCL version."""
		return self._runtime.get_device_cl_version(deviceID)

	def get_available_steppers(self) -> List[str]:
		"""Get the list of valid time stepper names."""
		return self._integrator.get_available_steppers()

	def get_program_string(self) -> str:
		"""Get the clODE OpenCL program string."""
		return self._integrator.get_program_string()

	def print_status(self) -> None:
		"""Print the simulator status info."""
		self._integrator.print_status()

	def print_devices(self) -> None:
		"""Print the available devices."""
		self._runtime.print_devices()


__all__ = ["Simulator", "SolverParams", "Stepper"]