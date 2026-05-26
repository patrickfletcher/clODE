from __future__ import annotations

from enum import Enum
from typing import Any, Callable, Dict, List, Mapping, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt

from .._opencl.executors import OpenCLTransientExecutor
from .._opencl.runtime import OpenCLRuntime
from ..problem._core import ProblemInfo
from ..problem.ivp import InitialValueProblem
from ..problem.python import OpenCLRhsEquation
from ..problem._core import RhsSource
from ..runtime import (
	CLDeviceType,
	CLVendor,
	OpenCLResource,
	_clode_root_dir,
	initialize_runtime,
)
from ..runtime.selection import RuntimeSelection
from ._state import SolverState, TransientCache
from .params import (
	_DEFAULT_INTEGRATION_SETTINGS,
	_DEFAULT_TRAJECTORY_OUTPUT_SETTINGS,
	SolverParams,
	_resolve_solver_params,
)


class Stepper(Enum):
	"""Supported time-stepping methods for the simulation backends."""

	euler = "euler"
	heun = "heun"
	rk4 = "rk4"
	bs23 = "bs23"
	dormand_prince = "dopri5"
	stochastic_euler = "seuler"


class Simulator:
	"""Base class for simulating an ensemble of instances of an ODE system.

	It provides the core functionality for advancing an ensemble in time without
	storing intermediate trajectory samples. Use it directly for final-state or
	transient-only workloads, or subclass it for feature extraction and trajectory
	storage workflows.
	"""

	_integrator: OpenCLTransientExecutor
	_runtime: OpenCLResource
	_runtime_selection: RuntimeSelection
	_single_precision: bool
	_stepper: Stepper
	_pi: ProblemInfo
	_rhs_source: RhsSource

	_cl_program_is_valid: bool = False

	_sp: SolverParams
	_solver_state: SolverState
	_transient_cache: TransientCache
	_ivp: InitialValueProblem

	_ensemble_size: int
	_ensemble_shape: Tuple

	def __init__(
		self,
		variables: Optional[Dict[str, float]] = None,
		parameters: Optional[Dict[str, float]] = None,
		aux: Optional[List[str]] = None,
		num_noise: int = 0,
		src_file: Optional[str] = None,
		rhs_equation: Optional[OpenCLRhsEquation] = None,
		supplementary_equations: List[Callable[[Any], Any]] | None = None,
		stepper: Stepper = Stepper.rk4,
		dt: float = _DEFAULT_INTEGRATION_SETTINGS.dt,
		dtmax: float = _DEFAULT_INTEGRATION_SETTINGS.dtmax,
		abstol: float = _DEFAULT_INTEGRATION_SETTINGS.abstol,
		reltol: float = _DEFAULT_INTEGRATION_SETTINGS.reltol,
		max_steps: int = _DEFAULT_INTEGRATION_SETTINGS.max_steps,
		max_store: int = _DEFAULT_TRAJECTORY_OUTPUT_SETTINGS.max_store,
		nout: int = _DEFAULT_TRAJECTORY_OUTPUT_SETTINGS.nout,
		solver_parameters: Optional[SolverParams] = None,
		t_span: Tuple[float, float] = (0.0, 1000.0),
		single_precision: bool = True,
		device_type: Optional[CLDeviceType] = None,
		vendor: Optional[CLVendor] = None,
		platform_id: Optional[int] = None,
		device_id: Optional[int] = None,
		ivp: Optional[InitialValueProblem] = None,
	) -> None:
		"""Create a simulator for one ODE model and one ensemble configuration.

		Args:
			variables: Mapping from state-variable name to its default initial value.
			parameters: Mapping from parameter name to its default value.
			aux: Ordered auxiliary-variable names written by the RHS.
			num_noise: Number of Wiener-process inputs expected by the RHS.
			src_file: Path to an OpenCL source file or XPP model file.
			rhs_equation: Typed Python RHS function to convert to OpenCL.
			supplementary_equations: Additional typed Python helper functions emitted
				into the generated OpenCL source before `rhs_equation`.
			stepper: Time-stepping method used by the backend.
			dt: Initial or fixed time step.
			dtmax: Maximum time step for adaptive steppers.
			abstol: Absolute tolerance for adaptive steppers.
			reltol: Relative tolerance for adaptive steppers.
			max_steps: Maximum number of integration steps per solve.
			max_store: Maximum number of stored time samples for trajectory solves.
			nout: Output stride for stored trajectories.
			solver_parameters: Optional prebuilt solver-parameter bundle. When
				provided, it overrides the scalar solver arguments above.
			t_span: Initial integration interval as `(t0, tf)`.
			single_precision: Whether to build the backend in single precision.
			device_type: Preferred OpenCL device class for runtime selection.
			vendor: Preferred OpenCL vendor for runtime selection.
			platform_id: Explicit OpenCL platform index.
			device_id: Explicit OpenCL device index on the selected platform.
			ivp: Optional explicit initial-value problem object. When provided, it
				replaces the model-definition arguments above.

		Raises:
			ValueError: If the model-definition arguments are inconsistent.
		"""

		self._ivp = self._coerce_initial_value_problem(
			ivp=ivp,
			variables=variables,
			parameters=parameters,
			aux=aux,
			num_noise=num_noise,
			src_file=src_file,
			rhs_equation=rhs_equation,
			supplementary_equations=supplementary_equations,
		)
		self._rhs_source = self._ivp.rhs_source
		self._pi = self._ivp.problem_info
		self._stepper = stepper
		self._single_precision = single_precision
		self._runtime = initialize_runtime(
			device_type,
			vendor,
			platform_id,
			device_id,
		)
		self._runtime_selection = self._runtime.runtime_selection

		self._create_integrator()
		self._build_cl_program()
		self._solver_state = SolverState()
		self._transient_cache = TransientCache()

		self._sp = _resolve_solver_params(
			solver_parameters=solver_parameters,
			dt=dt,
			dtmax=dtmax,
			abstol=abstol,
			reltol=reltol,
			max_steps=max_steps,
			max_store=max_store,
			nout=nout,
		)
		self.set_solver_parameters()

		self.set_tspan(t_span=t_span)
		self._sync_problem_data_from_ivp()

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
	def ivp(self) -> InitialValueProblem:
		"""The initial-value problem currently owned by the simulator."""
		return self._ivp

	def _create_opencl_runtime(self) -> OpenCLRuntime:
		return OpenCLRuntime.from_selection(self._runtime_selection)

	def _create_integrator(self) -> None:
		self._integrator = OpenCLTransientExecutor(
			self._pi,
			self._rhs_source,
			self._stepper.value,
			self._single_precision,
			self._create_opencl_runtime(),
			_clode_root_dir,
		)

	def _build_cl_program(self):
		self._integrator.build_cl()
		self._cl_program_is_valid = True

	def _ensure_cl_program(self) -> None:
		if not self._cl_program_is_valid:
			self._build_cl_program()

	def _invalidate_solution_cache(self) -> None:
		self._transient_cache.invalidate()
		self._solver_state.invalidate_results()

	def _invalidate_runtime_caches(self) -> None:
		self._invalidate_solution_cache()

	def _coerce_initial_value_problem(
		self,
		*,
		ivp: InitialValueProblem | None,
		variables: Dict[str, float] | None,
		parameters: Dict[str, float] | None,
		aux: List[str] | None,
		num_noise: int,
		src_file: str | None,
		rhs_equation: OpenCLRhsEquation | None,
		supplementary_equations: List[Callable[[Any], Any]] | None,
	) -> InitialValueProblem:
		if ivp is not None:
			if any(
				value is not None
				for value in (
					variables,
					parameters,
					aux,
					src_file,
					rhs_equation,
					supplementary_equations,
				)
			) or num_noise != 0:
				raise ValueError(
					"Cannot combine ivp with variables, parameters, or other model-definition arguments"
				)
			return ivp

		if variables is None or parameters is None:
			raise ValueError("Must specify either ivp or both variables and parameters")

		return InitialValueProblem(
			variables=variables,
			parameters=parameters,
			aux=aux,
			num_noise=num_noise,
			src_file=src_file,
			rhs_equation=rhs_equation,
			supplementary_equations=supplementary_equations,
		)

	def _sync_problem_data_from_ivp(self) -> None:
		self._ensemble_size = self._ivp.ensemble_size
		self._ensemble_shape = self._ivp.ensemble_shape
		self._set_problem_data(
			self._ivp.get_initial_state(),
			self._ivp.get_parameter_values(),
		)

	def _sync_ivp_from_device_problem_data(self) -> None:
		if not self._solver_state.problem_data_needs_pull:
			return
		initial_state = self._pull_initial_state_from_runtime()
		self._ivp._set_problem_data(
			initial_state,
			self._ivp.get_parameter_values(),
			ensemble_shape=self._ensemble_shape,
		)
		self._solver_state.mark_problem_data_synced()

	def _pull_initial_state_from_runtime(self) -> np.ndarray:
		initial_state = np.array(
			self._integrator.get_x0(), dtype=np.float64
		).reshape((self._ensemble_size, self.num_variables), order="F")
		self._ivp._set_initial_state(initial_state)
		self._solver_state.mark_problem_data_synced()
		return initial_state

	def set_repeat_ensemble(self, num_repeats: int) -> None:
		"""Create a 1D ensemble by repeating one parameter/state configuration.

		When the current ensemble size is 1, the current device state and parameter
		vector are broadcast to `num_repeats` instances. When resizing from a larger
		ensemble, the stored default values are used instead.

		Args:
			num_repeats: Number of independent copies to create.
		"""
		self._sync_ivp_from_device_problem_data()
		self._ivp.set_repeat_ensemble(num_repeats)
		self._sync_problem_data_from_ivp()

	def set_ensemble(
		self,
		variables: Optional[
			Union[np.ndarray, Mapping[str, Union[float, List[float], np.ndarray]]]
		] = None,
		parameters: Optional[
			Union[np.ndarray, Mapping[str, Union[float, List[float], np.ndarray]]]
		] = None,
	) -> None:
		"""Set or resize the ensemble by supplying state and/or parameter values.

		You may pass full `(ensemble_size, n)` arrays or dictionaries keyed by
		variable/parameter name. Dictionary values may be scalars or arrays with a
		shared shape; scalar values are broadcast across the new ensemble.

		Unspecified entries fall back to the stored defaults. When the ensemble size
		stays unchanged or expands from size 1, the current device initial state is
		reused before defaults are applied.

		Args:
			variables: Full initial-state array or mapping from variable name to
				scalar/array values.
			parameters: Full parameter array or mapping from parameter name to
				scalar/array values.

		Raises:
			ValueError: If both inputs are omitted, if names are unknown, or if the
				provided array shapes are incompatible.
		"""
		self._sync_ivp_from_device_problem_data()
		self._ivp.set_ensemble(variables=variables, parameters=parameters)
		self._sync_problem_data_from_ivp()

	def _set_problem_data(
		self, initial_state: np.ndarray, parameters: np.ndarray
	) -> None:
		"""Set both initial state and parameters at the same time."""
		self._ensemble_size = initial_state.shape[0]
		self._invalidate_runtime_caches()
		self._ivp._set_problem_data(
			initial_state,
			parameters,
			ensemble_shape=self._ensemble_shape,
		)
		self._integrator.set_problem_data(
			initial_state.flatten(order="F"),
			parameters.flatten(order="F"),
		)
		self._solver_state.reset_problem_time(self._ensemble_shape)

	def _set_parameters(self, parameters: np.ndarray) -> None:
		"""Set the ensemble parameters without changing ensemble size."""
		self._invalidate_runtime_caches()
		self._ivp._set_parameter_values(parameters)
		self._integrator.set_pars(parameters.flatten(order="F"))
		self._solver_state.reset_problem_time(self._ensemble_shape)

	def _set_initial_state(self, initial_state: np.ndarray) -> None:
		"""Set the initial state without changing ensemble size."""
		self._invalidate_runtime_caches()
		self._ivp._set_initial_state(initial_state)
		self._integrator.set_x0(initial_state.flatten(order="F"))
		self._solver_state.reset_problem_time(self._ensemble_shape)

	def set_tspan(self, t_span: tuple[float, float]) -> None:
		"""Set the integration interval used by subsequent solves.

		Args:
			t_span: Time interval as `(t0, tf)`.
		"""
		self._solver_state.set_requested_window(t_span)
		self._integrator.set_tspan(t_span)
		self._invalidate_runtime_caches()

	def get_tspan(self) -> tuple[float, float]:
		"""Return the integration interval currently stored on the device."""
		self._solver_state.set_requested_window(tuple(self._integrator.get_tspan()))
		return self._solver_state.t_span

	def shift_tspan(self) -> None:
		"""Advance the stored integration interval by one interval length."""
		self._integrator.shift_tspan()
		self._solver_state.set_requested_window(tuple(self._integrator.get_tspan()))
		self._invalidate_runtime_caches()

	def advance_tspan_to_attained_final_time(
		self,
		*,
		atol: float = 1e-12,
		rtol: float = 0.0,
	) -> tuple[float, float]:
		"""Advance the requested window so it starts at the attained final time.

		Unlike `shift_tspan()`, this preserves the current requested duration but
		starts the next window from the attained `tf` returned by the previous
		solve. This is the exact split-window continuation path when the ensemble
		shares one final time.

		Args:
			atol: Absolute tolerance used when checking that the ensemble agrees on
				one attained final time.
			rtol: Relative tolerance used when checking that the ensemble agrees on
				one attained final time.

		Returns:
			The new requested time window.

		Raises:
			ValueError: If no solve has been run yet or the ensemble reached
				different final times.
		"""
		self.get_final_time()
		next_tspan = self._solver_state.attained_final_time_window(
			atol=atol,
			rtol=rtol,
		)
		self.set_tspan(next_tspan)
		return next_tspan

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
		"""Update solver parameters and push them to the device.

		Args:
			solver_parameters: Optional complete parameter bundle. When provided, it
				replaces the current solver settings.
			dt: Initial or fixed time step.
			dtmax: Maximum time step for adaptive steppers.
			abstol: Absolute tolerance for adaptive steppers.
			reltol: Relative tolerance for adaptive steppers.
			max_steps: Maximum number of integration steps per solve.
			max_store: Maximum number of stored samples for trajectory solves.
			nout: Output stride for stored trajectories.
		"""
		if solver_parameters is not None:
			self._sp = solver_parameters.copy()
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
		self._invalidate_runtime_caches()

	def get_solver_parameters(self) -> SolverParams:
		"""Return the solver parameters currently stored on the device."""
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
		"""Run the simulator without storing trajectories or observer features.

		Args:
			t_span: Optional time interval override for this solve.
			update_x0: Whether to promote the final state to the next initial state
				after the solve.
			fetch_results: Whether to fetch and return the final state immediately.

		Returns:
			The final-state array when `fetch_results` is true, otherwise `None`.
		"""

		self._ensure_cl_program()

		if t_span is not None:
			self.set_tspan(t_span=t_span)

		self._integrator.transient()
		self._invalidate_solution_cache()

		if update_x0:
			self._integrator.shift_x0()
			self._solver_state.continue_problem_time()

		if fetch_results:
			return self.get_final_state()

	def get_initial_state(self) -> np.ndarray:
		"""Return the initial-state array with shape `(ensemble_size, num_variables)`."""
		if self._solver_state.problem_data_needs_pull:
			return self._pull_initial_state_from_runtime()
		return self._ivp.get_initial_state()

	def get_final_state(self) -> np.ndarray:
		"""Return the final-state array with shape `(ensemble_size, num_variables)`.

		Raises:
			ValueError: If no solve has been run yet.
		"""
		if self._transient_cache.final_state is None:
			final_state = self._integrator.get_xf()

			if final_state is None:
				raise ValueError("Must run a simulation before getting final state")

			self._transient_cache.final_state = np.array(final_state, dtype=np.float64).reshape(
				(self._ensemble_size, self.num_variables), order="F"
			)
		return self._transient_cache.final_state

	def get_dt(self) -> np.ndarray:
		"""Return per-instance continuation step sizes with shape matching the ensemble shape.

		For adaptive steppers this is the next step size the controller would attempt on
		a continued solve, not the width of the last accepted step. Use
		`get_last_accepted_dt()` for the latter.
		"""
		if self._solver_state.current_dt is None:
			self._solver_state.current_dt = np.array(
				self._integrator.get_dt(), dtype=np.float64
			).reshape(self._ensemble_shape, order="F")
		return self._solver_state.current_dt

	def get_status(self) -> np.ndarray:
		"""Return per-instance solve status codes with shape matching the ensemble shape.

		Values correspond to `clode.SolverStatus`.
		"""
		if self._solver_state.status is None:
			status = self._integrator.get_status()
			if len(status) == 0:
				raise ValueError("Must run a simulation before getting status")
			self._solver_state.status = np.array(
				status, dtype=np.int32
			).reshape(self._ensemble_shape, order="F")
		return self._solver_state.status

	def get_step_count(self) -> np.ndarray:
		"""Return per-instance accepted step counts with shape matching the ensemble shape."""
		if self._solver_state.step_count is None:
			step_count = self._integrator.get_step_count()
			if len(step_count) == 0:
				raise ValueError("Must run a simulation before getting step count")
			self._solver_state.step_count = np.array(
				step_count, dtype=np.uint64
			).reshape(self._ensemble_shape, order="F")
		return self._solver_state.step_count

	def get_last_accepted_dt(self) -> np.ndarray:
		"""Return the last accepted step width for each instance.

		Values have shape matching the ensemble shape. If a solve ends before any step
		is accepted, the reported value is zero.
		"""
		if self._solver_state.last_accepted_dt is None:
			accepted_dt = self._integrator.get_last_accepted_dt()
			if len(accepted_dt) == 0:
				raise ValueError("Must run a simulation before getting accepted step width")
			self._solver_state.last_accepted_dt = np.array(
				accepted_dt, dtype=np.float64
			).reshape(self._ensemble_shape, order="F")
		return self._solver_state.last_accepted_dt

	def get_final_time(self) -> np.ndarray:
		"""Return per-instance final times with shape matching the ensemble shape."""
		if self._solver_state.final_time is None:
			final_time = self._integrator.get_tf()
			if len(final_time) == 0:
				raise ValueError("Must run a simulation before getting final time")
			self._solver_state.final_time = np.array(
				final_time, dtype=np.float64
			).reshape(self._ensemble_shape, order="F")
			if (
				self._solver_state.problem_data_needs_pull
				and self._solver_state.current_time is None
			):
				self._solver_state.current_time = np.array(
					self._solver_state.final_time, dtype=np.float64, copy=True
				)
		return self._solver_state.final_time

	@property
	def platform_id(self) -> int:
		"""The concrete OpenCL platform index selected for this simulator."""
		return self._runtime.platform_id

	@property
	def device_id(self) -> int:
		"""The concrete OpenCL device index selected for this simulator."""
		return self._runtime.device_id

	@property
	def runtime_selection(self) -> RuntimeSelection:
		"""The concrete runtime selection bound to this simulator instance."""
		return self._runtime_selection

	@property
	def runtime_description(self) -> str:
		"""Human-readable description of the selected OpenCL runtime."""
		return self._runtime.describe()

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