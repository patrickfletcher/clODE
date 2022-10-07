import typing
import src.clode as _clode
import numpy as np
from enum import Enum


class Observer(Enum):
    basic = "basic"
    threshold_2 = "thresh2"


class ObserverOutput:

    def __init__(self, observer_params,
                 results_array: np.array,
                 num_result_features: int,
                 variables: list[str],
                 observer_type: Observer,
                 feature_names: list[str]):
        #print(type(observer_params))
        self._op = observer_params
        self._data = results_array
        self._num_result_features = num_result_features
        self._vars = variables
        self._observer_type = observer_type
        self._feature_names = feature_names

        # if self._observer_type == Observer.basic:
        #     shape = (self._num_result_features, len(results_array) // self._num_result_features)
        # # elif self._observer_type == Observer.threshold_2:
        # #     # Remove tbuffer, xbuffer and dxbuffer
        # #     start_idx = 3 + 3 * 2 * num_result_features
        # #     end_idx = start_idx +
        # #     results_array = results_array[]
        # else:
        #     raise NotImplementedError(f"{self._observer_type} not implemented!")

        shape = (self._num_result_features, len(results_array) // self._num_result_features)
        #print(f"Results array, {len(results_array)}, foo {shape}, bar {num_result_features}")
        self._data = results_array.reshape(shape).transpose()
        #print(self._data[0, :])

    # def get_var_max(self, var: str):
    #     if self._observer_type == Observer.basic:
    #         f_var_ix = self._op.f_var_ix
    #         observed_var = self._vars[f_var_ix]
    #         if var != observed_var:
    #             raise RuntimeException(f"Variable {var} not tracked "
    #                                    f"by observer {self._observer_type} "
    #                                    f"(only {observed_var} is tracked)")
    #         return self._data[:, 0:1]
    #     raise NotImplementedError(f"{self._observer_type} not implemented!")
    #
    # def get_var_min(self, var: str):
    #     if self._observer_type == Observer.basic:
    #         f_var_ix = self._op.f_var_ix
    #         observed_var = self._vars[f_var_ix]
    #         if var != observed_var:
    #             raise RuntimeException(f"Variable {var} not tracked "
    #                                    f"by observer {self._observer_type} "
    #                                    f"(only {observed_var} is tracked)")
    #         return self._data[:, 1:2]
    #     raise NotImplementedError(f"{self._observer_type} not implemented!")
    #
    # def get_var_avg(self, var: str):
    #     if self._observer_type == Observer.basic:
    #         f_var_ix = self._op.f_var_ix
    #         observed_var = self._vars[f_var_ix]
    #         if var != observed_var:
    #             raise RuntimeException(f"Variable {var} not tracked "
    #                                    f"by observer {self._observer_type} "
    #                                    f"(only {observed_var} is tracked)")
    #         return self._data[:, 2:3]
    #     # elif self._observer_type == Observer.threshold_2:
    #     #     offset = 3 +
    #     raise NotImplementedError(f"{self._observer_type} not implemented!")

    def _get_var(self, var: str, op: str):
        try:
            index = self._feature_names.index(f"{op} {var}")
            return self._data[:, index: index+1]
        except ValueError:
            raise NotImplementedError(f"{self._observer_type} does not track {op} {var}!")

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



