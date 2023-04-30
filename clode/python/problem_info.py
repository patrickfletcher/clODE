from .runtime import get_cpp

_clode = get_cpp()


class ProblemInfo:
    def __init__(
        self, src_file: str, vars: [str], pars: [str], aux: [str], num_noise: int
    ):
        _pi = _clode.problem_info(
            src_file, len(vars), len(pars), len(aux), num_noise, vars, pars, aux
        )
