from enum import Enum


class Stepper(Enum):
    euler = "euler"
    rk4 = "rk4"
    dormand_prince = "dopri5"
