

```python
import numpy as np
import matplotlib.pyplot as plt
import clode


def x_inf(v: float, vx: float, sx: float) -> float:
    return 1.0 / (1.0 + exp((vx - v) / sx))


def s_inf(c: float, k_s: float) -> float:
    c2: float = c * c
    return c2 / (c2 + k_s * k_s)


def lactotroph(
    t: float,
    x_: List[float],
    p_: List[float],
    dx_: List[float],
    aux_: List[float],
    w_: List[float],
) -> None:
    v: float = x_[0]
    n: float = x_[1]
    c: float = x_[2]

    gca: float = p_[0]
    gk: float = p_[1]
    gsk: float = p_[2]
    gleak: float = p_[3]
    cm: float = p_[4]
    e_leak: float = p_[5]
    tau_n: float = p_[6]
    k_c: float = p_[7]

    e_ca: float = 60
    e_k: float = -75

    vm: float = -20
    vn: float = -5
    sm: float = 12
    sn: float = 10

    f_c: float = 0.01
    alpha: float = 0.0015
    k_s: float = 0.4

    ica: float = gca * x_inf(v, vm, sm) * (v - e_ca)
    ik: float = gk * n * (v - e_k)
    isk: float = gsk * s_inf(c, k_s) * (v - e_k)
    ileak: float = gleak * (v - e_leak)
    current: float = ica + ik + isk + ileak

    dv: float = -current / cm

    dn: float = (x_inf(v, vn, sn) - n) / tau_n

    dc: float = -f_c * (alpha * ica + k_c * c)

    dx_[0] = dv
    dx_[1] = dn
    dx_[2] = dc


integrator = clode.FeatureSimulator(
    rhs_equation=lactotroph,
    supplementary_equations=[x_inf, s_inf],
    variables={
        "v": -60,
        "n": 0.1,
        "c": 0.1,
    },
    parameters={
        "gca": 1.5,
        "gk": 2.0,
        "gsk": 0.5,
        "gleak": 0.1,
        "cm": 10.0,
        "e_leak": -60,
        "tau_n": 15,
        "k_c": 0.1,
    },
    observer=clode.Observer.neighbourhood_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
)

# TODO magic that creates grid

parameters = []

integrator.set_ensemble(parameters=parameters)



```