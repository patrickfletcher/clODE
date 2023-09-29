from enum import Enum

import numpy as np

from .runtime import get_cpp

_clode = get_cpp()


class Observer(Enum):
    basic = "basic"
    basic_all_variables = "basicall"
    local_max = "localmax"
    neighbourhood_1 = "nhood1"
    neighbourhood_2 = "nhood2"
    threshold_2 = "thresh2"


class ObserverOutput:
    def __init__(
        self,
        observer_params,
        results_array: np.array,
        num_result_features: int,
        variables: list[str],
        observer_type: Observer,
        feature_names: list[str],
    ):
        # print(type(observer_params))
        self._op = observer_params
        self._data = results_array
        self._num_result_features = num_result_features
        self._vars = variables
        self._observer_type = observer_type
        self._feature_names = feature_names

        shape = (
            self._num_result_features,
            len(results_array) // self._num_result_features,
        )
        self._data = results_array.reshape(shape).transpose()

    def get_feature_names(self):
        return self._feature_names

    def _get_var(self, var: str, op: str):
        try:
            index = self._feature_names.index(f"{op} {var}")
            return self._data[:, index : index + 1]
        except ValueError:
            raise NotImplementedError(
                f"{self._observer_type} does not track {op} {var}!"
            )

    def get_var_max(self, var: str):
        return self._get_var(var, "max")

    def get_var_min(self, var: str):
        return self._get_var(var, "min")

    def get_var_mean(self, var: str):
        return self._get_var(var, "mean")

    def get_var_max_dt(self, var: str):
        return self.get_var_max(f"d{var}/dt")

    def get_var_min_dt(self, var: str):
        return self.get_var_min(f"d{var}/dt")

    def get_var_count(self, var: str):
        return self._get_var("count", var)
