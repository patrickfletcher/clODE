from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from .._opencl.executors import OpenCLFeatureExecutor
from ..observers.types import Observer, ObserverParams
from ..problem.ivp import InitialValueProblem
from ..problem.python import OpenCLRhsEquation
from ..runtime import CLDeviceType, CLVendor, _clode_root_dir
from ._state import FeatureCache
from .base import Simulator, Stepper
from .params import SolverParams
from .results import ObserverOutput


class FeatureSimulator(Simulator):
	"""Simulator that computes observer features and event data on the device."""

	_feature_cache: FeatureCache
	_integrator: OpenCLFeatureExecutor

	def __init__(
		self,
		variables: Optional[Dict[str, float]] = None,
		parameters: Optional[Dict[str, float]] = None,
		aux: Optional[List[str]] = None,
		num_noise: int = 0,
		src_file: Optional[str] = None,
		rhs_equation: Optional[OpenCLRhsEquation] = None,
		supplementary_equations: Optional[List[Callable[[Any], Any]]] = None,
		stepper: Stepper = Stepper.rk4,
		dt: float = 0.1,
		dtmax: float = 1.0,
		abstol: float = 1e-6,
		reltol: float = 1e-3,
		max_steps: int = 10000000,
		max_store: int = 10000000,
		nout: int = 1,
		solver_parameters: Optional[SolverParams] = None,
		t_span: Tuple[float, float] = (0.0, 1000.0),
		single_precision: bool = True,
		device_type: Optional[CLDeviceType] = None,
		vendor: Optional[CLVendor] = None,
		platform_id: Optional[int] = None,
		device_id: Optional[int] = None,
		ivp: Optional[InitialValueProblem] = None,
		observer: Observer = Observer.basic_all_variables,
		event_var: str = "",
		feature_var: str = "",
		observer_max_event_count: int = 100,
		observer_max_event_timestamps: int = 0,
		observer_min_x_amp: float = 0.0,
		observer_min_imi: float = 0.0,
		observer_neighbourhood_radius: float = 0.05,
		observer_x_up_thresh: float = 0.3,
		observer_x_down_thresh: float = 0.2,
		observer_dx_up_thresh: float = 0,
		observer_dx_down_thresh: float = 0,
		observer_eps_dx: float = 0.0,
		observer_parameters: Optional[ObserverParams] = None,
	) -> None:
		"""Create a feature-extraction simulator.

		This constructor accepts the same model-definition, solver, and runtime
		selection arguments as `Simulator`, plus configuration for the built-in
		observers used by `features()`.

		Args:
			observer: Built-in observer mode used during feature extraction.
			event_var: Variable name used for event detection when the observer
				requires one.
			feature_var: Variable name used for feature readout when the observer
				distinguishes detection and measurement variables.
			observer_max_event_count: Maximum number of events accumulated by the
				observer.
			observer_max_event_timestamps: Maximum number of event timestamps
				retained in the output.
			observer_min_x_amp: Minimum accepted event amplitude.
			observer_min_imi: Minimum inter-event interval.
			observer_neighbourhood_radius: Neighborhood radius for neighborhood-based
				observers.
			observer_x_up_thresh: Rising value threshold for threshold observers.
			observer_x_down_thresh: Falling value threshold for threshold observers.
			observer_dx_up_thresh: Rising derivative threshold for threshold
				observers.
			observer_dx_down_thresh: Falling derivative threshold for threshold
				observers.
			observer_eps_dx: Derivative tolerance used near threshold crossings.
			observer_parameters: Optional complete observer-parameter bundle. When
				provided, it overrides the individual `observer_*` arguments above.
		"""

		self._observer_type = observer
		self._feature_cache = FeatureCache()
		problem_variable_names = (
			ivp.variable_names if ivp is not None else list((variables or {}).keys())
		)

		event_var_idx = (
			problem_variable_names.index(event_var) if event_var != "" else 0
		)
		feature_var_idx = (
			problem_variable_names.index(feature_var) if feature_var != "" else 0
		)

		if observer_parameters is not None:
			self._op = observer_parameters
		else:
			self._op = ObserverParams(
				event_var_idx,
				feature_var_idx,
				observer_max_event_count,
				observer_max_event_timestamps,
				observer_min_x_amp,
				observer_min_imi,
				observer_neighbourhood_radius,
				observer_x_up_thresh,
				observer_x_down_thresh,
				observer_dx_up_thresh,
				observer_dx_down_thresh,
				observer_eps_dx,
			)

		super().__init__(
			variables=variables,
			parameters=parameters,
			src_file=src_file,
			rhs_equation=rhs_equation,
			supplementary_equations=supplementary_equations,
			aux=aux,
			num_noise=num_noise,
			t_span=t_span,
			stepper=stepper,
			single_precision=single_precision,
			dt=dt,
			dtmax=dtmax,
			abstol=abstol,
			reltol=reltol,
			max_steps=max_steps,
			max_store=max_store,
			nout=nout,
			solver_parameters=solver_parameters,
			device_type=device_type,
			vendor=vendor,
			platform_id=platform_id,
			device_id=device_id,
			ivp=ivp,
		)

	def _create_integrator(self) -> None:
		self._integrator = OpenCLFeatureExecutor(
			self._pi,
			self._rhs_source,
			self._stepper.value,
			self._observer_type.value,
			self._op,
			self._single_precision,
			self._create_opencl_runtime(),
			_clode_root_dir,
		)

	def _invalidate_feature_cache(self) -> None:
		self._feature_cache.invalidate()
		self._invalidate_solution_cache()

	def _invalidate_runtime_caches(self) -> None:
		super()._invalidate_runtime_caches()
		self._feature_cache.invalidate()

	def set_observer(self, observer_type: Observer) -> None:
		"""Switch to a different built-in observer.

		Args:
			observer_type: New observer mode to use for subsequent feature solves.
		"""
		if observer_type != self._observer_type:
			self._integrator.set_observer(observer_type.value)
			self._observer_type = observer_type
			self._cl_program_is_valid = False
			self._invalidate_feature_cache()

	def set_observer_parameters(
		self,
		op: Optional[ObserverParams] = None,
		event_var: Optional[str] = None,
		feature_var: Optional[str] = None,
		max_event_count: Optional[int] = None,
		max_event_timestamps: Optional[int] = None,
		min_amp: Optional[float] = None,
		min_imi: Optional[float] = None,
		nhood_radius: Optional[float] = None,
		x_up_threshold: Optional[float] = None,
		x_down_threshold: Optional[float] = None,
		dx_up_threshold: Optional[float] = None,
		dx_down_threshold: Optional[float] = None,
		eps_dx: Optional[float] = None,
	) -> None:
		"""Update observer parameters and push them to the device.

		Args:
			op: Optional complete observer-parameter bundle. When provided, it
				replaces the current observer settings.
			event_var: Variable name used for event detection.
			feature_var: Variable name used for feature readout.
			max_event_count: Maximum number of tracked events.
			max_event_timestamps: Maximum number of retained event timestamps.
			min_amp: Minimum accepted event amplitude.
			min_imi: Minimum inter-event interval.
			nhood_radius: Neighborhood radius for neighborhood-based observers.
			x_up_threshold: Rising value threshold for threshold observers.
			x_down_threshold: Falling value threshold for threshold observers.
			dx_up_threshold: Rising derivative threshold for threshold observers.
			dx_down_threshold: Falling derivative threshold for threshold observers.
			eps_dx: Derivative tolerance used near threshold crossings.
		"""
		current_max_event_timestamps = self._op.max_event_timestamps

		if op is not None:
			self._op = op
		else:
			if event_var is not None:
				self._op.e_var_ix = self.variable_names.index(event_var)
			if feature_var is not None:
				self._op.f_var_ix = self.variable_names.index(feature_var)
			if max_event_count is not None:
				self._op.max_event_count = max_event_count
			if max_event_timestamps is not None:
				self._op.max_event_timestamps = max_event_timestamps
			if min_amp is not None:
				self._op.min_amp = min_amp
			if min_imi is not None:
				self._op.min_imi = min_imi
			if nhood_radius is not None:
				self._op.nhood_radius = nhood_radius
			if x_up_threshold is not None:
				self._op.x_up_threshold = x_up_threshold
			if x_down_threshold is not None:
				self._op.x_down_threshold = x_down_threshold
			if dx_up_threshold is not None:
				self._op.dx_up_threshold = dx_up_threshold
			if dx_down_threshold is not None:
				self._op.dx_down_threshold = dx_down_threshold
			if eps_dx is not None:
				self._op.eps_dx = eps_dx

		if self._op.max_event_timestamps != current_max_event_timestamps:
			self._cl_program_is_valid = False

		self._integrator.set_observer_params(self._op)
		self._invalidate_feature_cache()

	def get_observer_parameters(self) -> ObserverParams:
		"""Return the current observer-parameter bundle from the backend."""
		return self._integrator.get_observer_params()

	def get_feature_names(self) -> List[str]:
		"""Get the list of feature names for the current observer."""
		return self._integrator.get_feature_names()

	def is_observer_initialized(self) -> bool:
		"""Return whether the active observer has completed its warmup pass."""
		return self._integrator.is_observer_initialized()

	def initialize_observer(self) -> None:
		"""Run the observer warmup pass, if the active observer requires one."""
		self._ensure_cl_program()
		self._integrator.initialize_observer()

	def features(
		self,
		t_span: Optional[Tuple[float, float]] = None,
		initialize_observer: Optional[bool] = None,
		update_x0: bool = True,
		fetch_results: bool = True,
	) -> Optional[ObserverOutput]:
		"""Run feature extraction for the current ensemble.

		Args:
			t_span: Optional time interval override for this solve.
			initialize_observer: Whether to rerun the observer initialization pass
				before extracting features. If `None`, use the backend default for the
				active observer.
			update_x0: Whether to promote the final state to the next initial state
				after the solve.
			fetch_results: Whether to fetch and return the observer output
				immediately.

		Returns:
			An `ObserverOutput` when `fetch_results` is true, otherwise `None`.
		"""
		self._ensure_cl_program()

		if t_span is not None:
			self.set_tspan(t_span=t_span)

		if initialize_observer is not None:
			self._integrator.features(initialize_observer)
		else:
			self._integrator.features()

		self._invalidate_solution_cache()
		self._feature_cache.mark_result_pending()

		if update_x0:
			self._integrator.shift_x0()
			self._solver_state.continue_problem_time()

		if fetch_results:
			return self.get_observer_results()

	def get_observer_results(self) -> ObserverOutput:
		"""Return the most recent observer output.

		Raises:
			ValueError: If `features()` has not been run yet.
		"""
		if not self._feature_cache.has_result:
			raise ValueError("Must run features() before getting observer results")

		if self._feature_cache.feature_array is None:
			self._feature_cache.feature_array = self._integrator.get_f()
			self._feature_cache.num_features = self._integrator.get_n_features()

		if (
			self._feature_cache.feature_array is None
			or self._feature_cache.num_features is None
		):
			raise ValueError("Must run features() before getting observer results")

		self._feature_cache.feature_array = np.array(
			self._feature_cache.feature_array, dtype=np.float64
		).reshape((self._ensemble_size, self._feature_cache.num_features), order="F")

		return ObserverOutput(
			self._op,
			self._feature_cache.feature_array,
			self._feature_cache.num_features,
			self.variable_names,
			self._observer_type,
			self._integrator.get_feature_names(),
			self._ensemble_shape,
		)


__all__ = ["FeatureSimulator"]