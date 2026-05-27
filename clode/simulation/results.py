from __future__ import annotations

from typing import Any, Optional, Tuple

import numpy as np
from numpy.lib import recfunctions as rfn

from ..observers.types import EventOutputSettings, Observer


class TrajectoryOutput:
    """Structured trajectory data returned by `TrajectorySimulator`.

    The `x`, `dx`, and `aux` arrays are exposed as structured NumPy arrays keyed by
    variable name. Use `to_ndarray()` when you need an unstructured numeric array.
    """

    def __init__(
        self,
        t: np.ndarray[Any, np.dtype[np.float64]],
        x: np.ndarray[Any, np.dtype[np.float64]],
        dx: np.ndarray[Any, np.dtype[np.float64]],
        aux: np.ndarray[Any, np.dtype[np.float64]],
        variable_names: list[str],
        aux_names: list[str],
    ) -> None:
        self.t = t

        x_dtype = np.dtype(
            {"names": variable_names, "formats": [np.float64] * len(variable_names)}
        )
        self.x = rfn.unstructured_to_structured(x, dtype=x_dtype)
        self.dx = rfn.unstructured_to_structured(dx, dtype=x_dtype)

        if len(aux_names) > 0:
            aux_dtype = np.dtype(
                {"names": aux_names, "formats": [np.float64] * len(aux_names)}
            )
            self.aux = rfn.unstructured_to_structured(aux, dtype=aux_dtype)

        self._variable_names = variable_names
        self._aux_names = aux_names

    def __repr__(self) -> str:
        return f"TrajectoryOutput( length: {len(self.t)}, variable names: {self._variable_names}, aux variable names: {self._aux_names} )"

    def to_ndarray(
        self, slot: str, **kwargs: Any
    ) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return one trajectory slot as a plain NumPy array.

        Args:
            slot: One of `"x"`, `"dx"`, or `"aux"`.

        Returns:
            An unstructured NumPy array for the requested slot.
        """

        if slot == "x":
            return rfn.structured_to_unstructured(self.x, **kwargs)
        if slot == "dx":
            return rfn.structured_to_unstructured(self.dx, **kwargs)
        if slot == "aux":
            return rfn.structured_to_unstructured(self.aux, **kwargs)
        raise ValueError(f"Unknown trajectory slot: {slot}")


class ObserverOutput:
    """Structured feature and event data returned by `FeatureSimulator`."""

    def __init__(
        self,
        event_output_settings: EventOutputSettings,
        feature_array: np.ndarray[Any, np.dtype[np.float64]],
        num_features: int,
        variables: list[str],
        observer_type: Observer,
        feature_names: list[str],
        ensemble_shape: Tuple,
    ) -> None:
        self._event_output_settings = event_output_settings
        self._num_features = num_features
        self._vars = variables
        self._observer_type = observer_type
        self._feature_names = feature_names
        self._ensemble_shape = ensemble_shape

        feature_dtype = np.dtype(
            {"names": feature_names, "formats": [np.float64] * len(feature_names)}
        )
        self.F = rfn.unstructured_to_structured(feature_array, dtype=feature_dtype)

    def __repr__(self) -> str:
        ensemble_size = len(self.F[self._feature_names[0]])
        num_features = len(self._feature_names)
        feature_names = self._feature_names
        return f"ObserverOutput( ensemble size: {ensemble_size}, number of features: {num_features}, feature_names: {feature_names})"

    def to_ndarray(self, **kwargs: Any) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return the feature table as an unstructured NumPy array."""

        return rfn.structured_to_unstructured(self.F, **kwargs)

    def get_feature_names(self) -> list[str]:
        """Return the ordered list of available feature names."""

        return self._feature_names

    def _get_var(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        try:
            result = self.F[var].squeeze().reshape(self._ensemble_shape)
            result = result[0] if result.size == 1 else result
            return result
        except ValueError as exc:
            raise NotImplementedError(f"{self._observer_type} does not track {var}!") from exc

    def get_var_max(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return the tracked maximum for one variable or derived quantity."""

        return self._get_var(" ".join(["max", var]))

    def get_var_min(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return the tracked minimum for one variable or derived quantity."""

        return self._get_var(" ".join(["min", var]))

    def get_var_mean(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return the tracked mean for one variable or derived quantity."""

        return self._get_var(" ".join(["mean", var]))

    def get_var_max_slope(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_max(f"d{var}/dt")

    def get_var_min_slope(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_min(f"d{var}/dt")

    def get_var_mean_slope(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        return self.get_var_mean(f"d{var}/dt")

    def get_var_count(self, var: str) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return the tracked count for one event or feature family."""

        return self._get_var(" ".join([var, "count"]))

    def get_event_data(
        self, name: str, type: Optional[str] = "time"
    ) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return stacked event data for one named event stream.

        Args:
            name: Event stream name such as `"up"`, `"down"`, or `"event"`.
            type: Event field to read, typically `"time"`.
        """

        event_features = [
            feature_name
            for feature_name in self._feature_names
            if name in feature_name
            and "event" in feature_name
            and "count" not in feature_name
        ]
        if len(event_features) == 0:
            raise NotImplementedError(
                f"{self._observer_type} does not track {name} event {type}s!"
            )
        data = []
        for event_idx in range(0, self._event_output_settings.max_event_timestamps):
            datapoint = self._get_var(f"{name} event {type} {event_idx}")
            if np.all(datapoint == 0):
                break
            data.append(datapoint)
        return np.stack(data, axis=-1).squeeze()

    def get_timestamps(
        self, var: str = "event"
    ) -> np.ndarray[Any, np.dtype[np.float64]]:
        """Return tracked event timestamps for one event family."""

        first_key = f"{var} event time 0"
        if first_key not in self._feature_names:
            raise NotImplementedError(
                f"{self._observer_type} does not track {var} event times!"
            )
        data = []
        for key_idx in range(0, self._event_output_settings.max_event_timestamps):
            datapoint = self._get_var(f"{var} event time {key_idx}")
            if np.all(datapoint == 0):
                break
            datapoint = (
                datapoint[np.newaxis] if len(datapoint.shape) == 0 else datapoint
            )
            data.append(datapoint)
        if data:
            return np.stack(data, axis=1).squeeze()
        return np.array([], dtype=np.float64)


__all__ = ["ObserverOutput", "TrajectoryOutput"]