from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from .._opencl.factory import create_trajectory_backend
from ..problem.python import OpenCLRhsEquation
from ..runtime import CLDeviceType, CLVendor, _clode_root_dir
from ._protocols import TrajectoryBackend
from .base import Simulator, Stepper
from .params import SolverParams
from .results import TrajectoryOutput


class TrajectorySimulator(Simulator):
	"""Simulator class that stores trajectories."""

	_device_t: np.ndarray[Any, np.dtype[np.float64]] | None
	_device_x: np.ndarray[Any, np.dtype[np.float64]] | None
	_device_dx: np.ndarray[Any, np.dtype[np.float64]] | None
	_device_aux: np.ndarray[Any, np.dtype[np.float64]] | None
	_integrator: TrajectoryBackend

	def __init__(
		self,
		variables: Dict[str, float],
		parameters: Dict[str, float],
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
		device_ids: Optional[List[int]] = None,
	) -> None:
		"""Construct a CLODE trajectory object."""

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
			device_ids=device_ids,
		)

		self._device_t = None
		self._device_x = None
		self._device_dx = None
		self._device_aux = None

	def _create_integrator(self) -> None:
		self._integrator = create_trajectory_backend(
			self._pi,
			self._rhs_source,
			self._stepper.value,
			self._single_precision,
			self._runtime,
			_clode_root_dir,
			runtime_selection=self._runtime_selection,
		)

	def trajectory(
		self,
		t_span: Optional[Tuple[float, float]] = None,
		update_x0: bool = True,
		fetch_results: bool = True,
	) -> Optional[List[TrajectoryOutput] | TrajectoryOutput]:
		"""Run a trajectory simulation."""
		self._ensure_cl_program()

		if t_span is not None:
			self.set_tspan(t_span=t_span)

		self._integrator.trajectory()
		self._device_t = self._device_x = self._device_dx = self._device_aux = None
		self._invalidate_solution_cache()

		if update_x0:
			self._integrator.shift_x0()
			self._device_initial_state = None

		if fetch_results:
			return self.get_trajectory()

	def get_trajectory(self) -> List[TrajectoryOutput] | TrajectoryOutput:
		"""Get the trajectory data."""

		self._device_n_stored = self._integrator.get_n_stored()
		self._device_t = self._integrator.get_t()
		self._device_x = self._integrator.get_x()
		self._device_dx = self._integrator.get_dx()
		self._device_aux = self._integrator.get_aux()

		if self._device_n_stored is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._device_t is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._device_x is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._device_dx is None:
			raise ValueError("Must run trajectory() before getting trajectory data")
		elif self._device_aux is None:
			raise ValueError("Must run trajectory() before getting trajectory data")

		t_shape = (self._ensemble_size, self._sp.max_store)
		self._device_t = np.array(
			self._device_t[: np.prod(t_shape)], dtype=np.float64
		).reshape(t_shape, order="F")

		data_shape = (self._ensemble_size, self.num_variables, self._sp.max_store)
		self._device_x = np.array(
			self._device_x[: np.prod(data_shape)], dtype=np.float64
		).reshape(data_shape, order="F")
		self._device_dx = np.array(
			self._device_dx[: np.prod(data_shape)], dtype=np.float64
		).reshape(data_shape, order="F")

		aux_shape = (self._ensemble_size, len(self.aux_names), self._sp.max_store)
		self._device_aux = np.array(
			self._device_aux[: np.prod(aux_shape)], dtype=np.float64
		).reshape(aux_shape, order="F")

		results = list()
		for i in range(self._ensemble_size):
			ni = self._device_n_stored[i] + 1
			ti = self._device_t[i, :ni].transpose()
			xi = self._device_x[i, :, :ni].transpose()
			dxi = self._device_dx[i, :, :ni].transpose()
			auxi = self._device_aux[i, :, :ni].transpose()
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