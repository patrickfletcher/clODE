from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np

from .function_converter import OpenCLRhsEquation
from .runtime import CLDeviceType, CLVendor, TrajectorySimulatorBase, _clode_root_dir
from .solver import Simulator, Stepper


class TrajectoryResult:
    t: np.ndarray[Any, np.dtype[np.float64]]
    x: np.ndarray[Any, np.dtype[np.float64]]

    def __init__(self, data: Dict[str, np.ndarray[Any, np.dtype[np.float64]]]) -> None:
        self.t = data["t"]
        self.x = data["x"]

    def __repr__(self) -> str:
        return f"TrajectoryResult(t={self.t}, x={self.x})"


class TrajectorySimulator(Simulator):
    _time_steps: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_time_steps: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_trajectories: np.ndarray[Any, np.dtype[np.float64]] | None
    _output_aux: np.ndarray[Any, np.dtype[np.float64]] | None
    _data: np.ndarray[Any, np.dtype[np.float64]] | None
    _aux: np.ndarray[Any, np.dtype[np.float64]] | None
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

        self._data = None
        self._output_trajectories = None
        self._time_steps = None
        self._output_time_steps = None
        self._output_aux = None
        self._aux = None

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
        self._output_time_steps = self._integrator.get_t()
        self._output_trajectories = self._integrator.get_x()

        # time_steps has one column per simulation (to support adaptive steppers)
        shape = (self._ensemble_size, self._max_store)
        arr = np.array(self._output_time_steps[: np.prod(shape)])
        self._time_steps = arr.reshape(shape, order="F").transpose((1, 0))

        data_shape = (self._ensemble_size, len(self.variable_names), self._max_store)
        arr = np.array(self._output_trajectories[: np.prod(data_shape)])
        self._data = arr.reshape(data_shape, order="F").transpose((2, 1, 0))

        # Check for None values
        if self._data is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._time_steps is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._output_time_steps is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._output_trajectories is None:
            raise ValueError("Must run trajectory() before getting trajectory data")
        elif self._n_stored is None:
            raise ValueError("Must run trajectory() before getting trajectory data")

        # list of trajectories, each stored as dict:
        results = list()
        for i in range(self._ensemble_size):
            ni = self._n_stored[i]
            ti = self._time_steps[:ni, i]
            xi = self._data[:ni, :, i]
            result = TrajectoryResult({"t": ti, "x": xi})
            results.append(result)

        return results

    def get_aux(self) -> List[np.ndarray[Any, np.dtype[np.float64]]]:
        """Get the auxiliary data.

        Returns:
            np.array: The auxiliary data.
        """
        _ = self.get_trajectory()
        self._output_aux = self._integrator.get_aux()

        shape = (self._ensemble_size, len(self.aux_variables), self._max_store)
        arr = np.array(self._output_aux[: np.prod(shape)])
        self._aux = arr.reshape(shape, order="F").transpose((2, 1, 0))

        results = list()
        for i in range(self._ensemble_size):
            ni = self._n_stored[i]
            aux_data = self._aux[:ni, :, i]
            results.append(aux_data)
        return results
