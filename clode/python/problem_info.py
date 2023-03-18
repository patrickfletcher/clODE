from .runtime import _get_clode

_clode = _get_clode()


class ProblemInfo:

    def __init__(self, src_file: str, vars: [str], pars: [str], aux: [str],
                 num_noise: int):
        _pi = _clode.problem_info(src_file, len(vars), len(pars), len(aux),
                                  num_noise, vars, pars, aux)
