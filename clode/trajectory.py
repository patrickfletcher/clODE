from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
# from numpy.typing import NDArray

from .function_converter import OpenCLRhsEquation
from .runtime import CLDeviceType, CLVendor, TrajectorySimulatorBase, _clode_root_dir
from .solver import Simulator, Stepper


class TrajectoryResult:
    def __init__(
            self, 
            t: np.ndarray[Any, np.dtype[np.float64]], 
            x: np.ndarray[Any, np.dtype[np.float64]], 
            dx: np.ndarray[Any, np.dtype[np.float64]], 
            aux: np.ndarray[Any, np.dtype[np.float64]], 
            variables: list[str], 
            aux_variables: list[str],
            ) -> None:
        self.t = t
        self._x = x
        self._dx = dx
        self._aux = aux
        self._variables = variables
        self._aux_variables = aux_variables

    def __repr__(self) -> str:
        return f"TrajectoryResult(length:{len(self.t)}, variables:{self._variables}, aux variables:{self._aux_variables})"

    # instead of the following, could simply use numpy structured arrays? similarly in features. numpy 1.24 supports 3.8+
    def x(self, var: str|None = None) -> np.ndarray[Any, np.dtype[np.float64]]:
        if var is None:
            return self._x
        
        try:
            index = self._variables.index(var)
            return self._x[:, index : index + 1].squeeze()
        except ValueError:
            raise NotImplementedError(
                f"{var} is not a valid variable name!"
            )
        
    def dx(self, var: str|None = None) -> np.ndarray[Any, np.dtype[np.float64]]:
        if var is None:
            return self._dx
        
        try:
            index = self._variables.index(var)
            return self._dx[:, index : index + 1].squeeze()
        except ValueError:
            raise NotImplementedError(
                f"{var} is not a valid variable name!"
            )

    def aux(self, var: str|None = None) -> np.ndarray[Any, np.dtype[np.float64]]:
        if var is None:
            return self._aux
        
        try:
            index = self._aux_variables.index(var)
            return self._aux[:, index : index + 1].squeeze()
        except ValueError:
            raise NotImplementedError(
                f"{var} is not a valid aux variable name!"
            )
        
class TrajectorySimulator(Simulator):
    _time_steps: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_t: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_x: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_dx: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_aux: np.ndarray[Any, np.dtype[np.float64]] | None
    _x_data: np.ndarray[Any, np.dtype[np.float64]] | None
    _aux_data: np.ndarray[Any, np.dtype[np.float64]] | None
    _integrator: TrajectorySimulatorBase

    def __init__(
        self,
        variables: Dict[str, float],
        parameters: Dict[str, float],
        src_file: Optional[str] = None,
        rhs_equation: Optional[OpenCLRhsEquation] = None,
        supplementary_equations: Optional[List[Callable[[Any], Any]]] = None,
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        t_span: Tuple[float, float] = (0.0, 1000.0),
        stepper: Stepper = Stepper.rk4,
        single_precision: bool = True,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-4,
        max_steps: int = 1000000,
        max_store: int = 1000000,
        nout: int = 1,
        device_type: CLDeviceType | None = None,
        vendor: CLVendor | None = None,
        platform_id: int | None = None,
        device_id: int | None = None,
        device_ids: List[int] | None = None,
    ) -> None:
        # We use Google-style docstrings: https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html
        """Initialize a CLODE trajectory object.

        Args:
            src_file (str): The path to the source file to be simulated.  If the file ends with ".xpp", it will be converted to a CLODE source file.
            variable_names (List[str]): The names of the variables to be simulated.
            parameter_names (List[str]): The names of the parameters to be simulated.
            aux (Optional[List[str]], optional): The names of the auxiliary variables to be simulated. Defaults to None.
            num_noise (int, optional): The number of noise variables to be simulated. Defaults to 0.
            tspan (Tuple[float, float], optional): The time span to simulate over. Defaults to (0.0, 1000.0).
            stepper (Stepper, optional): The stepper to use. Defaults to Stepper.rk4.
            single_precision (bool, optional): Whether to use single precision. Defaults to True.
            dt (float, optional): The initial time step. Defaults to 0.1.
            dtmax (float, optional): The maximum time step. Defaults to 1.0.
            abstol (float, optional): The absolute tolerance. Defaults to 1e-6.
            reltol (float, optional): The relative tolerance. Defaults to 1e-3.
            max_steps (int, optional): The maximum number of steps. Defaults to 1000000.
            max_store (int, optional): The maximum number of time steps to store. Defaults to 1000000.
            nout (int, optional): The number of output time steps. Defaults to 1.
            device_type (Optional[CLDeviceType], optional): The type of device to use. Defaults to None.
            vendor (Optional[CLVendor], optional): The vendor of the device to use. Defaults to None.
            platform_id (Optional[int], optional): The platform ID of the device to use. Defaults to None.
            device_id (Optional[int], optional): The device ID of the device to use. Defaults to None.
            device_ids (Optional[List[int]], optional): The device IDs of the devices to use. Defaults to None.

        Raises:
            ValueError: If the source file does not exist.

        Returns (CLODETrajectory): The initialized CLODE trajectory object.
        """

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
            device_type=device_type,
            vendor=vendor,
            platform_id=platform_id,
            device_id=device_id,
            device_ids=device_ids,
        )

        self._x_data = None
        self._time_steps = None
        self._output_t = None
        self._output_x = None
        self._output_dx = None
        self._output_aux = None
        self._aux_data = None

    def _build_integrator(self) -> None:
        self._integrator = TrajectorySimulatorBase(
            self._pi,
            self._stepper.value,
            self._single_precision,
            self._runtime,
            _clode_root_dir,
        )

    def trajectory(self, update_x0: bool = True) -> List[TrajectoryResult]:
        """Run a trajectory simulation.

        Returns:
            List[TrajectoryResult]
        """
        if not self.is_initialized:
            raise RuntimeError("Simulator is not initialized")

        self._integrator.trajectory()
        if update_x0:
            self.shift_x0()

        return self.get_trajectory()

    def get_trajectory(self) -> List[TrajectoryResult]:
        """Get the trajectory data.

        Returns:
            np.array: The trajectory data.
        """

        # fetch data from device
        self._n_stored = self._integrator.get_n_stored()
        self._output_t = self._integrator.get_t()
        self._output_x = self._integrator.get_x()
        self._output_dx = self._integrator.get_dx()
        self._output_aux = self._integrator.get_aux()

        # Check for None values
        if self._n_stored is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._output_t is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._output_x is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._output_dx is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._output_aux is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        
        # time_steps has one column per simulation (to support adaptive steppers)
        shape = (self._ensemble_size, self._max_store)
        arr = np.array(self._output_t[: np.prod(shape)])
        self._time_steps = arr.reshape(shape, order="F").transpose((1, 0))

        data_shape = (self._ensemble_size, len(self.variable_names), self._max_store)
        arr = np.array(self._output_x[: np.prod(data_shape)])
        self._x_data = arr.reshape(data_shape, order="F").transpose((2, 1, 0))
        arr = np.array(self._output_dx[: np.prod(data_shape)])
        self._dx_data = arr.reshape(data_shape, order="F").transpose((2, 1, 0))

        aux_shape = (self._ensemble_size, len(self.aux_variables), self._max_store)
        arr = np.array(self._output_aux[: np.prod(aux_shape)])
        self._aux_data = arr.reshape(aux_shape, order="F").transpose((2, 1, 0))

        # list of trajectories, each stored as dict:
        results = list()
        for i in range(self._ensemble_size):
            ni = self._n_stored[i]
            ti = self._time_steps[:ni, i]
            xi = self._x_data[:ni, :, i]
            dxi = self._dx_data[:ni, :, i]
            auxi = self._aux_data[:ni, :, i]
            result = TrajectoryResult(t=ti, x=xi, dx=dxi, aux=auxi, variables=self.variable_names, aux_variables=self.aux_variables)
            results.append(result)

        return results