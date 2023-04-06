# Feature extraction

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

### XPP

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
    return x**(1 / 3.)


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

    parameters = [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] + \
                 list(range(5, end))

    x0 = np.tile([1, 1], (len(parameters), 1))

    pars_v = np.array([[par] for par in parameters])
    integrator.initialize(x0, pars_v)

    integrator.transient()
    integrator.features()
    observer_output = integrator.get_observer_results()
    
    return observer_output

vdp_dormand_prince(100, "vdp_oscillator.xpp")
```
