from __future__ import annotations

from typing import List, Optional, Tuple

import numpy as np

from .runtime import _clode_root_dir, get_cpp
from .solver import CLODE, Stepper

_clode = get_cpp()


# TrajectoryOutput?
# We have a collection of num_simulations trajectories, (t, X), stacked into matrices. Each trajectory may have a different number of total stored time steps.
# - for convenience, would be nice to have access patterns something like:
# >>> trajectory_output.t[0] --> t[: nstored[0], 0]
# >>> trajectory_output.X[0] --> X[: nstored[0], :, 0])
# >>> trajectory_output.X[var, 0] --> X[: nstored[0], var, 0]
# >>> trajectory_output.t, .X --> (t[: max(nstored), :], X[: max(nstored), :, :])
# class TrajectoryOutput:


class CLODETrajectory(CLODE):
    def __init__(
        self,
        src_file: str,
        variable_names: List[str],
        parameter_names: List[str],
        aux: Optional[List[str]] = None,
        num_noise: int = 0,
        tspan: Tuple[float, float] = (0.0, 1000.0),
        stepper: Stepper = Stepper.rk4,
        single_precision: bool = True,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-4,
        max_steps: int = 1000000,
        max_store: int = 1000000,
        nout: int = 1,
        device_type: _clode.cl_device_type | None = None,
        vendor: _clode.cl_vendor | None = None,
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
            device_type (Optional[_clode.cl_device_type], optional): The type of device to use. Defaults to None.
            vendor (Optional[_clode.cl_vendor], optional): The vendor of the device to use. Defaults to None.
            platform_id (Optional[int], optional): The platform ID of the device to use. Defaults to None.
            device_id (Optional[int], optional): The device ID of the device to use. Defaults to None.
            device_ids (Optional[List[int]], optional): The device IDs of the devices to use. Defaults to None.

        Raises:
            ValueError: If the source file does not exist.

        Returns (CLODETrajectory): The initialized CLODE trajectory object.
        """

        super().__init__(
            src_file=src_file,
            variable_names=variable_names,
            parameter_names=parameter_names,
            aux=aux,
            num_noise=num_noise,
            tspan=tspan,
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

        self._integrator = _clode.clode_trajectory(
            self._pi, stepper.value, single_precision, self._runtime, _clode_root_dir
        )
        self._integrator.build_cl()

    def initialize(
        self,
        x0: np.array,
        parameters: np.array,
        tspan: Tuple[float, float] | None = None,
        seed: int | None = None,
    ) -> None:
        """Initialize the trajectory object.

        Args:
            x0 (np.array): The initial conditions.
            parameters (np.array): The parameters.
            tspan (Tuple[float, float], optional): The time span to simulate over. Defaults to None.
            seed (int, optional): The seed for the random number generator. Defaults to None.

        Raises:
            ValueError: If the initial conditions or parameters are not the correct shape.

        Returns:
            None
        """

        if len(x0.shape) != 2:
            raise ValueError("Must provide rows of initial variables")

        if x0.shape[1] != len(self.vars):
            raise ValueError(
                f"Length of initial condition vector {len(x0.shape[1])}"
                f" does not match number of variables {len(self.vars)}"
            )

        if len(parameters.shape) != 2:
            raise ValueError("Must provide rows of parameters")

        if parameters.shape[1] != len(self.pars):
            raise ValueError(
                f"Length of parameters vector {parameters.shape[1]}"
                f" does not match number of parameters {len(self.pars)}"
            )

        self._data = None
        self._time_steps = None
        self._number_of_simulations = parameters.shape[0]

        if tspan is not None:
            self.tspan = tspan

        self._integrator.initialize(
            self.tspan,
            x0.transpose().flatten(),
            parameters.transpose().flatten(),
            self._sp,
        )
        self.seed_rng(seed)

    def trajectory(self) -> None:
        """Run a trajectory simulation.

        Returns:
            None
        """
        self._integrator.trajectory()

    def get_trajectory(self) -> np.array:
        """Get the trajectory data.

        Returns:
            np.array: The trajectory data.
        """

        # fetch data from device
        self._n_stored = self._integrator.get_n_stored()
        self._output_time_steps = self._integrator.get_t()
        self._output_trajectories = self._integrator.get_x()

        # time_steps has one column per simulation (to support adaptive steppers)
        shape = (self._number_of_simulations, self._max_store)
        arr = np.array(self._output_time_steps[: np.prod(shape)])
        self._time_steps = arr.reshape(shape, order="F").transpose((1, 0))

        shape = (self._number_of_simulations, len(self.vars), self._max_store)
        arr = np.array(self._output_trajectories[: np.prod(shape)])
        self._data = arr.reshape(shape, order="F").transpose((2, 1, 0))

        # list of trajectories, each stored as dict:
        result = list()
        for i in range(self._number_of_simulations):
            ni = self._n_stored[i]
            ti = self._time_steps[:ni, i]
            xi = self._data[:ni, :, i]
            result.append({"t": ti, "X": xi})

        return result
