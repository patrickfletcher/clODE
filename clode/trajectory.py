from __future__ import annotations

from typing import Any, Callable, Dict, List, Optional, Tuple

import numpy as np
# from numpy.typing import NDArray
from numpy.lib import recfunctions as rfn

from .function_converter import OpenCLRhsEquation
from .runtime import CLDeviceType, CLVendor, _clode_root_dir
from clode.cpp.clode_cpp_wrapper import TrajectorySimulatorBase
from .solver import Simulator, Stepper


class TrajectoryOutput:
    def __init__(
            self, 
            t: np.ndarray[Any, np.dtype[np.float64]], 
            x: np.ndarray[Any, np.dtype[np.float64]], 
            dx: np.ndarray[Any, np.dtype[np.float64]], 
            aux: np.ndarray[Any, np.dtype[np.float64]], 
            variables: list[str], 
            aux_names: list[str],
            ) -> None:
        
        self.t = t

        x_dtype = np.dtype({"names":variables, "formats":[np.float64]*len(variables)})
        self.x = rfn.unstructured_to_structured(x, dtype=x_dtype)
        self.dx = rfn.unstructured_to_structured(dx, dtype=x_dtype)

        if len(aux_names)>0:
            aux_dtype = np.dtype({"names":aux_names, "formats":[np.float64]*len(aux_names)})
            self.aux = rfn.unstructured_to_structured(aux, dtype=aux_dtype)

        self._variables = variables
        self._aux_names = aux_names

    def __repr__(self) -> str:
        return f"TrajectoryOutput( length: {len(self.t)}, variables: {self._variables}, aux variables: {self._aux_names} )"
    
    # helper to convert back to unstructured ndarray
    # --> make this a class property decorator?
    # alternatively: self.x.view(np.float64).reshape(-1,len(variables))?
    def to_ndarray(self, slot:str, **kwargs):
        if slot=="x":
            return rfn.structured_to_unstructured(self.x, **kwargs)
        elif slot=="dx":
            return rfn.structured_to_unstructured(self.dx, **kwargs)
        elif slot=="aux":
            return rfn.structured_to_unstructured(self.aux, **kwargs)
        
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
            t_span (Tuple[float, float], optional): The time span to simulate over. Defaults to (0.0, 1000.0).
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

        self._max_store = self._sp.max_store
        self._n_out = self._sp.nout

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
        self._integrator.build_cl()
        self._cl_program_is_valid = True

    # TODO[feature]: chunk time - keep max_store to a reasonable level (device-dependent), loop solve/get until t_span is covered.
    def trajectory(self, 
                   t_span:Optional[Tuple[float, float]] = None, 
                   update_x0: bool = True, 
                   fetch_results:bool = True,
        ) -> Optional[List[TrajectoryOutput]|TrajectoryOutput]:
        """Run a trajectory simulation.
        
        Args:
        t_span (tuple[float, float]): Time interval for integration.
        update_x0 (bool): After the simulation, whether to overwrite the initial state buffer with the final state
        fetch_results (bool): Whether to fetch the feature results from the device and return them here

        Returns:
            List[TrajectoryOutput]
        """

        if not self.is_initialized:
            raise RuntimeError("Simulator is not initialized")

        if t_span is not None:
            self.set_tspan(t_span=t_span)

        self._integrator.trajectory()
        # invalidate _device_t, _device_x, _device_dx, _device_aux, _device_final_state
        # self._device_t = self._device_x = self._device_dx = self._device_aux = self._device_final_state = None

        if update_x0:
            self._integrator.shift_x0()
            # invalidate _device_initial_state

        if fetch_results:
            return self.get_trajectory()

    # TODO: specialize? support individual getters too
    def get_trajectory(self) -> List[TrajectoryOutput]|TrajectoryOutput:
        """Get the trajectory data.

        Returns:
            TrajectoryOutput
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
        t_shape = (self._ensemble_size, self._max_store)
        arr = np.array(self._output_t[: np.prod(t_shape)])
        self._time_steps = arr.reshape(t_shape, order="F").transpose() 

        data_shape = (self._ensemble_size, len(self.variable_names), self._max_store)
        arr = np.array(self._output_x[: np.prod(data_shape)])
        self._x_data = arr.reshape(data_shape, order="F").transpose((2, 1, 0))
        arr = np.array(self._output_dx[: np.prod(data_shape)])
        self._dx_data = arr.reshape(data_shape, order="F").transpose((2, 1, 0))

        aux_shape = (self._ensemble_size, len(self.aux_names), self._max_store)
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
            result = TrajectoryOutput(t=ti, x=xi, dx=dxi, aux=auxi, variables=self.variable_names, aux_names=self.aux_names)
            results.append(result)

        return results[0] if self._ensemble_size==1 else results