import runtime

_clode = runtime.get_cpp()


class SolverParams:
    def __init__(
        self,
        dt: float = 0.1,
        dtmax: float = 1.0,
        abstol: float = 1e-6,
        reltol: float = 1e-3,
        max_steps: int = 1000000,
        max_store: int = 1000000,
        nout: int = 1,
    ):
        self._sp = _clode.solver_params(
            dt, dtmax, abstol, reltol, max_steps, max_store, nout
        )
