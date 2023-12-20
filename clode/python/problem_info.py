from __future__ import annotations

from typing import List

from .runtime import get_cpp

_clode = get_cpp()


class ProblemInfo:
    def __init__(
        self,
        src_file: str,
        vars: List[str],
        pars: List[str],
        aux: List[str],
        num_noise: int,
    ):
        self._pi = _clode.problem_info(
            src_file, len(vars), len(pars), len(aux), num_noise, vars, pars, aux
        )
