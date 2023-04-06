# CLODE - Getting started

## Installation

See [installation](installation.md) for instructions on how to install CLODE.

## Usage (Feature extraction)

CLODE can simulate a large number of ODEs simultaneously using OpenCL.

You will need an OpenCL Right-Hand-Side (RHS) function.

This is an implementation of tne Van der Pol oscillator:

```c++
void getRHS(const realtype t,
            const realtype var[],
            const realtype par[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {
    realtype m = par[0];
    realtype w = par[1];
    realtype k = par[2];
    realtype H = par[3];

    realtype y1 = var[0];
    realtype y2 = var[1];

    realtype dy1 = y2;
    realtype dy2 = (w - k * y2) / m;

    derivatives[0] = dy1;
    derivatives[1] = dy2;
    aux[0] = y1 - H;
}
```

Note that you must specify a function called getRHS, with the above prototype.

You can declare other valid OpenCL functions in the same file, and
use them inside getRHS.

In Python this is straightforward:

```python
import clode

tspan = (0.0, 1000.0)

# Create a feature extractor
integrator = clode.CLODEFeatures(
    # This is your source file. 
    src_file="van_der_pol_oscillator.cl",
    # The variables must match getRHS
    variable_names=["x", "y"],
    parameter_names=["mu"],
    num_noise=1, # Number of Weiner noise variables (at least one)
    observer=clode.Observer.threshold_2, # Choose an observer
    stepper=clode.Stepper.dormand_prince, # Choose a stepper
    tspan=tspan,
)

# Set the parameters
# In this example, they are repeated for each simulation,
# but in real-world useage, they will be different
parameters = [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0] +
             list(range(5, end))

x0 = np.tile([1, 1], (len(parameters), 1))

pars_v = np.array([[par] for par in parameters])
integrator.initialize(x0, pars_v)

# Run the simulation for tspan time, 
# to get rid of transient effects
integrator.transient()

# Run the simulation for tspan time,
# to get the features
integrator.features()

# Get the results from the feature observer
observer_output = integrator.get_observer_results()

```

Note: There is a bug in the noise generation. It is currently not possible to
set the number of noise variables to 0. Set to 1 or more.

# Implementation details

The Python library wraps a CPP library, clode_cpp_wrapper.[so|dll]
The CPP library assumes that the variables/parameters are grouped
by columns, i.e. if your variables are a, b and c,
the CPP library expects data in the format [aaaabbbbcccc].
The Python library expects data in the format
[[a, b, c], [a, b, c], [a, b, c], ...].