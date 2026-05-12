from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(slots=True)
class ProblemInfo:
    """Static shape information for a modeled ODE problem.

    Attributes:
        src_file: Human-readable source label for the problem definition.
        vars: Ordered state-variable names.
        pars: Ordered parameter names.
        aux: Ordered auxiliary-variable names.
        num_noise: Number of Wiener-process inputs.
    """

    src_file: str = ""
    vars: list[str] = field(default_factory=list)
    pars: list[str] = field(default_factory=list)
    aux: list[str] = field(default_factory=list)
    num_noise: int = 1

    def __post_init__(self) -> None:
        self.vars = list(self.vars)
        self.pars = list(self.pars)
        self.aux = list(self.aux)
        self.num_noise = int(self.num_noise)
        if self.num_noise < 0:
            raise ValueError("num_noise must be non-negative")

    @property
    def num_var(self) -> int:
        return len(self.vars)

    @property
    def num_par(self) -> int:
        return len(self.pars)

    @property
    def num_aux(self) -> int:
        return len(self.aux)


__all__ = ["ProblemInfo"]