from math import exp


class OpenCLFunction:
    def to_string(self) -> str:
        raise NotImplementedError


class OpenCLExp(OpenCLFunction):
    @staticmethod
    def __call__(arg1: float) -> float:
        return exp(arg1)

    def __str__(self) -> str:
        return "exp"


class OpenCLMin(OpenCLFunction):
    @staticmethod
    def __call__(arg1: float, arg2: float) -> float:
        return min(arg1, arg2)

    def __str__(self) -> str:
        return "fmin"
