# Feature extraction

...Under Construction...

CLODE can simulate a large number of ODEs simultaneously using OpenCL.
To keep memory usage low, CLODE extracts features on-the-fly
using configurable observers. These observers are written in OpenCL
and can capture properties like the minimum, maximum, or average
of a variable.

Advanced observers can also capture local maxima, neighbourhoods,
and thresholds.

## Observers

CLODE's observers are highly configurable. You can choose the following:

* basic - Captures one variable
* basic_all_variables - Captures all variables
* local_max - Captures local maxima
* neighbourhood_2
* threshold_2 - Captures all values above a threshold

## Example

The following example extracts features from the Van der Pol oscillator
using the dormand_prince45 integrator.

```python
import numpy as np
import matplotlib.pyplot as plt
import clode


# Van der Pol Dormand Prince oscillator
def getRHS(
        t: float,
        var: list[float],
        par: list[float],
        derivatives: list[float],
        aux: list[float],
        wiener: list[float],
) -> None:
    mu: float = par[0]
    x: float = var[0]
    y: float = var[1]

    dx: float = y
    dy: float = mu * (1 - x ** 2) * y - x

    derivatives[0] = dx
    derivatives[1] = dy


def scipy_solve_ivp_wrapper(func, aux=None, wiener=None):
    if aux is None:
        aux = []
    if wiener is None:
        wiener = []

    def wrapper(t, y, *args):
        dydt = np.zeros_like(y)
        func(t, y, args, dydt, aux, wiener)
        return dydt

    return wrapper


wrap = scipy_solve_ivp_wrapper(getRHS)
xx = solve_ivp(
    wrap,
    [0, 1000],
    [1, 1],
    args=(-1, 0, 1),
    atol=1e-10,
    rtol=1e-10,
    mxstep=1000000,
)

# Invoke the wrapper with scipy's odeint

getRHS = scipy_odeint_wrapper(getRHS)

from scipy.integrate import odeint

res = odeint(
    getRHS,
    [1, 1],
    [0, 1000],
    args=([-1], [], []),
    atol=1e-10,
    rtol=1e-10,
    mxstep=1000000,
)

integrator = clode.FeatureSimulator(
    rhs_equation=getRHS,
    variable_names=["x", "y"],
    parameter_names=["mu"],
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.dormand_prince,
    tspan=(0.0, 1000.0),
)

parameters = [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

x0 = np.tile([1, 1], (len(parameters), 1))

pars_v = np.array([[par] for par in parameters])
integrator.set_ensemble(x0, pars_v)

integrator.transient()
integrator.features()
observer_output = integrator.get_observer_results()

periods = observer_output.get_var_max("period")
plt.plot(parameters, periods[:, 0])
plt.title("Van der Pol oscillator")
plt.xlabel("mu")
plt.ylabel("period")

plt.show()
```

## New definition

```python
import numpy as np
import matplotlib.pyplot as plt
import clode


# Van der Pol Dormand Prince oscillator
def get_rhs(
        t: float,
        var: list[float],
        mu: float,
        kl: float,
        weiner1: float,
) -> list[float]:
    x: float = var[0]
    y: float = var[1]

    kk: float = weiner1 * kl

    dx: float = y
    dy: float = mu * (1 - x ** 2) * y - x

    aux: float = x + kk

    return [dx, dy]


ivp = clode.IVP(
    rhs=get_rhs,
    variables: dict[str, float] = {"x": 1.0, "y": 1.0},
parameters: dict[str, float] = {"mu": 0.1},
aux: list["str"] = ["aux"],
noise: list["str"] = ["weiner1"]
)


integrator = clode.FeatureSimulator(
    ivp=ivp,
    solver=clode.Solver.dormand_prince,
)

integrator.set_ensemble(
    parameters={"mu": [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]},
    variables={"x": [1.0, 2.0], },
)

integrator.transient()
integrator.features()
observer_output = integrator.get_observer_results()

periods = observer_output.get_var_max("period")

plt.plot(parameters, periods[:, 0])
plt.title("Van der Pol oscillator")
plt.xlabel("mu")
plt.ylabel("period")

plt.show()
```

## XPP

```xpp
init x = 0.1 y = 0.1

par mu = 0.1
y' = mu * (1 - x*x) * y - x
x' = y

@ dt=0.05, total=5000, maxstor=20000000
@ bounds=10000000, xp=t, yp=v
@ xlo=0, xhi=5000, ylo=-75, yhi=0
@ method=Euler

```

### Python

```python
import numpy as np

import clode


def cuberoot(x):
    return x ** (1 / 3.)


def vdp_dormand_prince(end: int, input_file: str):
    tspan = (0.0, 1000.0)

    integrator = clode.CLODEFeatures(
        src_file=input_file,
        variable_names=["x", "y"],
        parameter_names=["mu"],
        num_noise=0,
        observer=clode.Observer.threshold_2,
        stepper=clode.Stepper.dormand_prince,
        tspan=tspan,
    )

    parameters = [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] +
    list(range(5, end))


x0 = np.tile([1, 1], (len(parameters), 1))

pars_v = np.array([[par] for par in parameters])
integrator.set_ensemble(x0, pars_v)

integrator.transient()
integrator.features()
observer_output = integrator.get_observer_results()

return observer_output

vdp_dormand_prince(100, "vdp_oscillator.xpp")
```

```python
def scipy_solve_ivp_wrapper(func, aux=None, wiener=None):
    if aux is None:
        aux = []
    if wiener is None:
        wiener = []
    def wrapper(t, y, *args):
        dydt = np.zeros_like(y)
        func(t, y, args, dydt, aux, wiener)
        return dydt

    return wrapper

wrap = scipy_solve_ivp_wrapper(getRHS)
xx = solve_ivp(
    wrap,
    [0, 1000],
    [1, 1],
    args=(-1, 0, 1),
    atol=1e-10,
    rtol=1e-10,
    mxstep=1000000,
)

# Invoke the wrapper with scipy's odeint

getRHS = scipy_odeint_wrapper(getRHS)
```
