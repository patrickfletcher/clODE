from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any, Callable

import numpy as np

from .definition import ProblemInfo
from .python import OpenCLConverter, OpenCLRhsEquation
from .source import RhsSource, create_rhs_source, load_rhs_source
from .xpp import convert_xpp_file

ArrayValue = float | Sequence[float] | np.ndarray[Any, np.dtype[np.float64]]


def _prepare_rhs_source(
	src_file: str | None = None,
	rhs_equation: OpenCLRhsEquation | None = None,
	supplementary_equations: Sequence[Callable[[Any], Any]] | None = None,
) -> tuple[RhsSource, OpenCLRhsEquation | None]:
	if src_file is not None and rhs_equation is not None:
		raise ValueError("Cannot specify both src_file and rhs_equation")
	if src_file is not None:
		input_file = convert_xpp_file(src_file) if src_file.endswith(".xpp") else src_file
		return load_rhs_source(input_file), None
	if rhs_equation is not None:
		converter = OpenCLConverter()
		if supplementary_equations is not None:
			for equation in supplementary_equations:
				converter.convert_to_opencl(equation)
		eqn = converter.convert_to_opencl(
			rhs_equation, mutable_args=[3, 4], function_name="getRHS"
		)
		return create_rhs_source("clode_rhs.cl", eqn), rhs_equation
	raise ValueError("Must specify either src_file or rhs_equation")


class InitialValueProblem:
	"""Semantic owner of one ODE problem definition and its batch-shaped inputs."""

	_problem_info: ProblemInfo
	_rhs_source: RhsSource
	_python_rhs_equation: OpenCLRhsEquation | None
	_variable_defaults: dict[str, float]
	_parameter_defaults: dict[str, float]
	_initial_state: np.ndarray[Any, np.dtype[np.float64]]
	_parameter_values: np.ndarray[Any, np.dtype[np.float64]]
	_ensemble_size: int
	_ensemble_shape: tuple[int, ...]

	def __init__(
		self,
		variables: Mapping[str, float],
		parameters: Mapping[str, float],
		aux: Sequence[str] | None = None,
		num_noise: int = 0,
		src_file: str | None = None,
		rhs_equation: OpenCLRhsEquation | None = None,
		supplementary_equations: Sequence[Callable[[Any], Any]] | None = None,
	) -> None:
		if aux is None:
			aux = []

		self._variable_defaults = dict(variables)
		self._parameter_defaults = dict(parameters)
		self._rhs_source, self._python_rhs_equation = _prepare_rhs_source(
			src_file=src_file,
			rhs_equation=rhs_equation,
			supplementary_equations=supplementary_equations,
		)
		self._problem_info = ProblemInfo(
			self._rhs_source.origin_label,
			list(self._variable_defaults.keys()),
			list(self._parameter_defaults.keys()),
			list(aux),
			num_noise,
		)

		default_initial_state = np.array(
			list(self._variable_defaults.values()), dtype=np.float64, ndmin=2
		)
		default_parameters = np.array(
			list(self._parameter_defaults.values()), dtype=np.float64, ndmin=2
		)
		self._set_problem_data(default_initial_state, default_parameters, ensemble_shape=(1,))

	@property
	def problem_info(self) -> ProblemInfo:
		"""Derived static shape metadata used by the runtime implementation."""
		return self._problem_info

	@property
	def rhs_source(self) -> RhsSource:
		"""Prepared RHS source bundle used by the OpenCL build path."""
		return self._rhs_source

	@property
	def variable_names(self) -> list[str]:
		return list(self._problem_info.vars)

	@property
	def parameter_names(self) -> list[str]:
		return list(self._problem_info.pars)

	@property
	def aux_names(self) -> list[str]:
		return list(self._problem_info.aux)

	@property
	def num_variables(self) -> int:
		return self._problem_info.num_var

	@property
	def num_parameters(self) -> int:
		return self._problem_info.num_par

	@property
	def num_aux(self) -> int:
		return self._problem_info.num_aux

	@property
	def num_noise(self) -> int:
		return self._problem_info.num_noise

	@property
	def ensemble_size(self) -> int:
		return self._ensemble_size

	@property
	def ensemble_shape(self) -> tuple[int, ...]:
		return self._ensemble_shape

	@property
	def default_initial_state(self) -> np.ndarray[Any, np.dtype[np.float64]]:
		return np.array(list(self._variable_defaults.values()), dtype=np.float64)

	@property
	def default_parameters(self) -> np.ndarray[Any, np.dtype[np.float64]]:
		return np.array(list(self._parameter_defaults.values()), dtype=np.float64)

	@property
	def python_rhs_equation(self) -> OpenCLRhsEquation | None:
		"""Original Python RHS callable when this IVP was authored in Python."""
		return self._python_rhs_equation

	def get_initial_state(self) -> np.ndarray[Any, np.dtype[np.float64]]:
		"""Return the current batch-shaped initial-state array."""
		return np.array(self._initial_state, dtype=np.float64, copy=True)

	def get_parameter_values(self) -> np.ndarray[Any, np.dtype[np.float64]]:
		"""Return the current batch-shaped parameter array."""
		return np.array(self._parameter_values, dtype=np.float64, copy=True)

	def _set_problem_data(
		self,
		initial_state: np.ndarray[Any, np.dtype[np.float64]],
		parameters: np.ndarray[Any, np.dtype[np.float64]],
		*,
		ensemble_shape: tuple[int, ...] | None = None,
	) -> None:
		initial_state_array = np.array(initial_state, dtype=np.float64, copy=True, ndmin=2)
		parameter_array = np.array(parameters, dtype=np.float64, copy=True, ndmin=2)

		if initial_state_array.ndim != 2 or initial_state_array.shape[1] != self.num_variables:
			raise ValueError(
				f"initial_state must be a matrix with {self.num_variables} columns"
			)
		if parameter_array.ndim != 2 or parameter_array.shape[1] != self.num_parameters:
			raise ValueError(
				f"parameters must be a matrix with {self.num_parameters} columns"
			)
		if initial_state_array.shape[0] != parameter_array.shape[0]:
			raise ValueError("initial_state and parameters must have the same size")

		new_size = int(initial_state_array.shape[0])
		if ensemble_shape is None:
			if new_size == self._ensemble_size:
				ensemble_shape = self._ensemble_shape
			else:
				ensemble_shape = (new_size,)
		if int(np.prod(ensemble_shape)) != new_size:
			raise ValueError("ensemble_shape must match the ensemble size")

		self._initial_state = initial_state_array
		self._parameter_values = parameter_array
		self._ensemble_size = new_size
		self._ensemble_shape = ensemble_shape

	def _set_initial_state(
		self, initial_state: np.ndarray[Any, np.dtype[np.float64]]
	) -> None:
		self._set_problem_data(
			initial_state,
			self._parameter_values,
			ensemble_shape=self._ensemble_shape,
		)

	def _set_parameter_values(
		self, parameters: np.ndarray[Any, np.dtype[np.float64]]
	) -> None:
		self._set_problem_data(
			self._initial_state,
			parameters,
			ensemble_shape=self._ensemble_shape,
		)

	def set_repeat_ensemble(self, num_repeats: int) -> None:
		"""Repeat one state/parameter configuration into a 1D ensemble."""
		initial_state, parameters = self._make_problem_data(
			new_size=num_repeats,
			new_shape=(num_repeats, 1),
		)
		self._set_problem_data(initial_state, parameters, ensemble_shape=(num_repeats, 1))

	def set_ensemble(
		self,
		variables: np.ndarray[Any, np.dtype[np.float64]]
		| Mapping[str, ArrayValue]
		| None = None,
		parameters: np.ndarray[Any, np.dtype[np.float64]]
		| Mapping[str, ArrayValue]
		| None = None,
	) -> None:
		"""Set or resize the IVP batch using state and parameter values."""
		if variables is None and parameters is None:
			raise ValueError("initial_state and parameters cannot both be None")

		validated_variables = self._validate_ensemble_input(
			values=variables,
			names=self.variable_names,
			label="variables",
			expected_width=self.num_variables,
		)
		validated_parameters = self._validate_ensemble_input(
			values=parameters,
			names=self.parameter_names,
			label="parameters",
			expected_width=self.num_parameters,
		)

		var_size = 1
		var_shape = (1,)
		if isinstance(validated_variables, np.ndarray):
			var_size = validated_variables.shape[0]
			var_shape = (var_size, 1)
		elif isinstance(validated_variables, Mapping):
			var_size, var_shape = self._shape_from_mapping(validated_variables, "variables")

		par_size = 1
		par_shape = (1,)
		if isinstance(validated_parameters, np.ndarray):
			par_size = validated_parameters.shape[0]
			par_shape = (par_size, 1)
		elif isinstance(validated_parameters, Mapping):
			par_size, par_shape = self._shape_from_mapping(validated_parameters, "parameters")

		if var_size > 1 and par_size > 1 and var_size != par_size:
			raise ValueError(
				"Arrays specified for parameters and initial states must have the same size"
			)

		new_size = var_size if var_size > 1 else par_size
		new_shape = var_shape if var_size > 1 else par_shape
		initial_state, parameter_values = self._make_problem_data(
			variables=validated_variables,
			parameters=validated_parameters,
			new_size=new_size,
			new_shape=new_shape,
		)
		self._set_problem_data(initial_state, parameter_values, ensemble_shape=new_shape)

	def evaluate_rhs(
		self,
		t: float,
		y: Sequence[float] | np.ndarray[Any, np.dtype[np.float64]],
		*args: float,
		parameters: Mapping[str, float]
		| Sequence[float]
		| np.ndarray[Any, np.dtype[np.float64]]
		| None = None,
		aux: Sequence[float] | np.ndarray[Any, np.dtype[np.float64]] | None = None,
		wiener: Sequence[float] | np.ndarray[Any, np.dtype[np.float64]] | None = None,
	) -> np.ndarray[Any, np.dtype[np.float64]]:
		"""Evaluate a Python-authored RHS as a SciPy-style derivative function."""
		if self._python_rhs_equation is None:
			raise ValueError(
				"SciPy-style RHS evaluation is only available for IVPs built from a Python RHS"
			)

		state = np.asarray(y, dtype=np.float64)
		if state.ndim != 1 or state.size != self.num_variables:
			raise ValueError(
				f"y must be a one-dimensional array with {self.num_variables} entries"
			)

		parameter_values = self._resolve_parameter_values(args=args, parameters=parameters)
		aux_values = self._resolve_vector(aux, self.num_aux, "aux")
		wiener_values = self._resolve_vector(wiener, self.num_noise, "wiener")
		derivatives = np.zeros(self.num_variables, dtype=np.float64)

		self._python_rhs_equation(
			float(t),
			state,
			parameter_values,
			derivatives,
			aux_values,
			wiener_values,
		)
		return derivatives

	def __call__(
		self,
		t: float,
		y: Sequence[float] | np.ndarray[Any, np.dtype[np.float64]],
		*args: float,
	) -> np.ndarray[Any, np.dtype[np.float64]]:
		"""Delegate to `evaluate_rhs` so the IVP can be used with `solve_ivp`."""
		return self.evaluate_rhs(t, y, *args)

	def _shape_from_mapping(
		self,
		values: Mapping[str, np.ndarray[Any, np.dtype[np.float64]]],
		label: str,
	) -> tuple[int, tuple[int, ...]]:
		sizes = [value.size for value in values.values() if value.size > 1]
		shapes = [value.shape for value in values.values() if value.size > 1]
		if len(set(shapes)) > 1:
			shape_map = {key: value.shape for key, value in values.items() if value.size > 1}
			raise ValueError(f"Shape of arrays for {label} don't match: {shape_map}")
		if sizes:
			return sizes[0], shapes[0]
		return 1, (1,)

	def _make_problem_data(
		self,
		variables: Mapping[str, np.ndarray[Any, np.dtype[np.float64]]] | np.ndarray[Any, np.dtype[np.float64]] | None = None,
		parameters: Mapping[str, np.ndarray[Any, np.dtype[np.float64]]] | np.ndarray[Any, np.dtype[np.float64]] | None = None,
		new_size: int | None = None,
		new_shape: tuple[int, ...] | None = None,
	) -> tuple[np.ndarray[Any, np.dtype[np.float64]], np.ndarray[Any, np.dtype[np.float64]]]:
		if new_size is None or new_shape is None:
			raise ValueError("new_size and new_shape are required")
		if len(new_shape) == 1:
			new_shape = (new_size, 1)

		previous_size = self._ensemble_size
		valid_previous_size = previous_size == new_size or previous_size == 1

		if valid_previous_size:
			initial_state_array = self.get_initial_state()
			parameter_array = self.get_parameter_values()
		else:
			initial_state_array = np.array(self.default_initial_state, dtype=np.float64, ndmin=2)
			parameter_array = np.array(self.default_parameters, dtype=np.float64, ndmin=2)

		if initial_state_array.shape[0] == 1:
			initial_state_array = np.tile(initial_state_array, (new_size, 1))
		if parameter_array.shape[0] == 1:
			parameter_array = np.tile(parameter_array, (new_size, 1))

		if isinstance(variables, np.ndarray):
			initial_state_array = variables
		elif isinstance(variables, Mapping):
			for key, value in variables.items():
				index = self.variable_names.index(key)
				value_array = np.repeat(value, new_size) if value.size == 1 else value
				initial_state_array[:, index] = np.asarray(value_array).reshape(-1)

		if isinstance(parameters, np.ndarray):
			parameter_array = parameters
		elif isinstance(parameters, Mapping):
			for key, value in parameters.items():
				index = self.parameter_names.index(key)
				value_array = np.repeat(value, new_size) if value.size == 1 else value
				parameter_array[:, index] = np.asarray(value_array).reshape(-1)

		self._ensemble_size = new_size
		self._ensemble_shape = new_shape
		return initial_state_array, parameter_array

	def _resolve_parameter_values(
		self,
		*,
		args: tuple[float, ...],
		parameters: Mapping[str, float]
		| Sequence[float]
		| np.ndarray[Any, np.dtype[np.float64]]
		| None,
	) -> np.ndarray[Any, np.dtype[np.float64]]:
		if args and parameters is not None:
			raise ValueError("Use either positional args or parameters=, not both")
		if args:
			if len(args) != self.num_parameters:
				raise ValueError(
					f"Expected {self.num_parameters} positional parameter overrides, got {len(args)}"
				)
			return np.asarray(args, dtype=np.float64)
		if parameters is None:
			return self.default_parameters
		if isinstance(parameters, Mapping):
			unknown = set(parameters.keys()) - set(self.parameter_names)
			if unknown:
				raise ValueError(f"Unknown parameter name(s): {unknown}")
			parameter_values = self.default_parameters
			for key, value in parameters.items():
				index = self.parameter_names.index(key)
				parameter_values[index] = float(value)
			return parameter_values
		parameter_values = np.asarray(parameters, dtype=np.float64)
		if parameter_values.ndim != 1 or parameter_values.size != self.num_parameters:
			raise ValueError(
				f"parameters must be a one-dimensional array with {self.num_parameters} entries"
			)
		return parameter_values

	def _resolve_vector(
		self,
		values: Sequence[float] | np.ndarray[Any, np.dtype[np.float64]] | None,
		expected_size: int,
		label: str,
	) -> np.ndarray[Any, np.dtype[np.float64]]:
		if values is None:
			return np.zeros(expected_size, dtype=np.float64)
		array = np.asarray(values, dtype=np.float64)
		if array.ndim != 1 or array.size != expected_size:
			raise ValueError(
				f"{label} must be a one-dimensional array with {expected_size} entries"
			)
		return array

	def _validate_ensemble_input(
		self,
		*,
		values: np.ndarray[Any, np.dtype[np.float64]] | Mapping[str, ArrayValue] | None,
		names: list[str],
		label: str,
		expected_width: int,
	) -> np.ndarray[Any, np.dtype[np.float64]] | Mapping[str, np.ndarray[Any, np.dtype[np.float64]]] | None:
		if isinstance(values, np.ndarray):
			if values.ndim != 2 or values.shape[1] != expected_width:
				raise ValueError(
					f"{label} must be a matrix with {expected_width} columns"
				)
			return np.array(values, dtype=np.float64, copy=True)
		if isinstance(values, Mapping):
			unknown = set(values.keys()) - set(names)
			if unknown:
				raise ValueError(f"Unknown {label[:-1]} name(s): {unknown}")
			return {
				key: np.asarray(value, dtype=np.float64)
				for key, value in values.items()
			}
		if values is not None:
			raise ValueError(
				f"Expected np.ndarray or Mapping for {label}, but got {type(values)}"
			)
		return None


__all__ = ["InitialValueProblem"]