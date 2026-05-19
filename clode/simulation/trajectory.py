from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from .._opencl.executors import OpenCLTrajectoryExecutor
from ..problem.ivp import InitialValueProblem
from ..problem.python import OpenCLRhsEquation
from ..runtime import CLDeviceType, CLVendor, _clode_root_dir
from ._state import TrajectoryCache
from .base import Simulator, Stepper
from .params import SolverParams, _TrajectoryOutputSettings
from .results import TrajectoryOutput


class TrajectorySimulator(Simulator):
	"""Simulator that stores time samples and returns `TrajectoryOutput` objects."""

	_trajectory_cache: TrajectoryCache
	_trajectory_output_settings: _TrajectoryOutputSettings
	_integrator: OpenCLTrajectoryExecutor

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
		reltol: float = 1e-4,
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
		ivp: Optional[InitialValueProblem] = None,
	) -> None:
		"""Create a trajectory simulator.

		This constructor accepts the same model-definition, solver, and runtime
		selection arguments as `Simulator`. In trajectory workflows, `max_store`
		and `nout` determine how many time samples are retained and how densely they
		are stored.
		"""
		initial_solver_parameters = (
			solver_parameters.copy()
			if solver_parameters is not None
			else SolverParams(dt, dtmax, abstol, reltol, max_steps, max_store, nout)
		)
		self._trajectory_cache = TrajectoryCache()
		self._trajectory_output_settings = (
			initial_solver_parameters.trajectory_output_settings
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
		self._integrator = OpenCLTrajectoryExecutor(
			self._pi,
			self._rhs_source,
			self._stepper.value,
			self._single_precision,
			self._create_opencl_runtime(),
			_clode_root_dir,
		)

	def _invalidate_runtime_caches(self) -> None:
		super()._invalidate_runtime_caches()
		self._trajectory_cache.invalidate()

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
		"""Update integration and trajectory-output settings.

		Integration changes invalidate solver state. Output-only changes only
		invalidate stored trajectory results.
		"""
		if all(
			value is None
			for value in (
				solver_parameters,
				dt,
				dtmax,
				abstol,
				reltol,
				max_steps,
				max_store,
				nout,
			)
		):
			super().set_solver_parameters()
			self._trajectory_output_settings = self._sp.trajectory_output_settings
			return

		previous_solver_parameters = self._sp.copy()
		previous_integration = previous_solver_parameters.integration_settings
		previous_output = self._trajectory_output_settings

		if solver_parameters is not None:
			updated_solver_parameters = solver_parameters.copy()
		else:
			updated_solver_parameters = previous_solver_parameters.copy()
			if dt is not None:
				updated_solver_parameters.dt = dt
			if dtmax is not None:
				updated_solver_parameters.dtmax = dtmax
			if abstol is not None:
				updated_solver_parameters.abstol = abstol
			if reltol is not None:
				updated_solver_parameters.reltol = reltol
			if max_steps is not None:
				updated_solver_parameters.max_steps = max_steps
			if max_store is not None:
				updated_solver_parameters.max_store = max_store
			if nout is not None:
				updated_solver_parameters.nout = nout

		if updated_solver_parameters == previous_solver_parameters:
			return

		updated_integration = updated_solver_parameters.integration_settings
		updated_output = updated_solver_parameters.trajectory_output_settings

		self._sp = updated_solver_parameters
		self._trajectory_output_settings = updated_output
		self._integrator.set_solver_params(self._sp)

		if updated_integration != previous_integration:
			self._invalidate_runtime_caches()
		elif updated_output != previous_output:
			self._trajectory_cache.invalidate()

	def trajectory(
		self,
		t_span: Optional[Tuple[float, float]] = None,
		update_x0: bool = True,
		fetch_results: bool = True,
	) -> Optional[List[TrajectoryOutput] | TrajectoryOutput]:
		"""Run a solve that stores trajectory samples.

		Args:
			t_span: Optional time interval override for this solve.
			update_x0: Whether to promote the final state to the next initial state
				after the solve.
			fetch_results: Whether to fetch and return the stored trajectory output
				immediately.

		Returns:
			A single `TrajectoryOutput` for a size-1 ensemble, a list of
			`TrajectoryOutput` objects for larger ensembles, or `None` when
			`fetch_results` is false.
		"""
		self._ensure_cl_program()

		if t_span is not None:
			self.set_tspan(t_span=t_span)

		self._integrator.trajectory()
		self._invalidate_solution_cache()
		self._trajectory_cache.mark_result_pending()

		if update_x0:
			self._integrator.shift_x0()
			self._solver_state.continue_problem_time()

		if fetch_results:
			return self.get_trajectory()

	def get_trajectory(self) -> List[TrajectoryOutput] | TrajectoryOutput:
		"""Return the most recently stored trajectory output.

		Returns:
			A single `TrajectoryOutput` for a size-1 ensemble or a list of outputs
			for larger ensembles.

		Raises:
			ValueError: If `trajectory()` has not been run yet.
		"""
		if not self._trajectory_cache.has_result:
			raise ValueError("Must run trajectory() before getting trajectory data")

		self._trajectory_cache.n_stored = np.asarray(
			self._integrator.get_n_stored(), dtype=np.int32
		)
		self._trajectory_cache.t = np.asarray(self._integrator.get_t(), dtype=np.float64)
		self._trajectory_cache.x = np.asarray(self._integrator.get_x(), dtype=np.float64)
		self._trajectory_cache.dx = np.asarray(self._integrator.get_dx(), dtype=np.float64)
		self._trajectory_cache.aux = np.asarray(self._integrator.get_aux(), dtype=np.float64)

		if self._trajectory_cache.n_stored is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._trajectory_cache.t is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._trajectory_cache.x is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._trajectory_cache.dx is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._trajectory_cache.aux is None:
			raise ValueError("Must run trajectory() before getting trajectory data")

		t_shape = (
			self._ensemble_size,
			self._trajectory_output_settings.max_store,
		)
		self._trajectory_cache.t = np.array(
			self._trajectory_cache.t[: np.prod(t_shape)], dtype=np.float64
		).reshape(t_shape, order="F")

		data_shape = (
			self._ensemble_size,
			self.num_variables,
			self._trajectory_output_settings.max_store,
		)
		self._trajectory_cache.x = np.array(
			self._trajectory_cache.x[: np.prod(data_shape)], dtype=np.float64
		).reshape(data_shape, order="F")
		self._trajectory_cache.dx = np.array(
			self._trajectory_cache.dx[: np.prod(data_shape)], dtype=np.float64
		).reshape(data_shape, order="F")

		aux_shape = (
			self._ensemble_size,
			len(self.aux_names),
			self._trajectory_output_settings.max_store,
		)
		self._trajectory_cache.aux = np.array(
			self._trajectory_cache.aux[: np.prod(aux_shape)], dtype=np.float64
		).reshape(aux_shape, order="F")

		results = list()
		for i in range(self._ensemble_size):
			ni = self._trajectory_cache.n_stored[i] + 1
			ti = self._trajectory_cache.t[i, :ni].transpose()
			xi = self._trajectory_cache.x[i, :, :ni].transpose()
			dxi = self._trajectory_cache.dx[i, :, :ni].transpose()
			auxi = self._trajectory_cache.aux[i, :, :ni].transpose()
			result = TrajectoryOutput(
				t=ti,
				x=xi,
				dx=dxi,
				aux=auxi,
				variable_names=self.variable_names,
				aux_names=self.aux_names,
			)
			results.append(result)

		return results[0] if self._ensemble_size == 1 else results


__all__ = ["TrajectorySimulator"]