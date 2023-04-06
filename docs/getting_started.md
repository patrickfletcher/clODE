# CLODE - Getting started

## Installation

See [installation](install.md) for instructions on how to install CLODE.

Once installed, you can have clODE query your system and display information about available OpenCL devices - see [Querying system OpenCL capabilities](querying_opencl.md)

## Usage Example - computing the period of the Van der Pol oscillator

[The Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator) can be written as a system of two differential equations:

$
\dot{x} = y\\
\dot{y} = \mu(1-x^2)y - x
$

Oscillations occur when $\mu>0$. Suppose we wish to measure the period of oscillations as $\mu$ varies. First, we will need to implement the vector field above as function - the right-hand-side (RHS) function - with the signature expected by clODE:

```c
void getRHS(const realtype t,
            const realtype x_[],
            const realtype p_[],
            realtype dx_[],
            realtype aux_[],
            const realtype w_[]) {

    /* State variables */
    realtype x = x_[0];
    realtype y = x_[1];

    /* Parameters */
    realtype mu = p_[0];

    /* Differential equations */
    realtype dx = y;
    realtype dy = mu * (1 - x*x) * y - x;

    /* Differential outputs */
    dx_[0] = dx;
    dx_[1] = dy;
}
```

Note that this is a simple C-language function. The ```realtype``` type declaration is a macro that will expand to ```float``` or ```double```, depending on our choice of single or double precision, which defaults to single precision.  Save this function in a file called ```van_der_pol_oscillator.cl```

This function supports additional use cases not used in this example: ```t``` for time-dependent ODE terms (non-autonomous systems), ```aux``` for auxiliary readout variables, and stochastic terms via ```w```, which provides Wiener variables ($w \sim Normal \; (0,dt)$). Finally, one can also declare other simple C functions in the same file before ```getRHS``, and use them inside getRHS. For further example use cases, see [TODO]

Next we will use a python script to define our parameters and set up the numerical simulation. Here we use clODE's feature detection mode - several features of the ODE solution, including the period of oscillation, will be measured "on the fly", without storing the trajectory itself.

```python
import clode
import numpy as np

# time span for our simulation
tspan = (0.0, 1000.0)

# Create the clODE feature extractor
integrator = clode.CLODEFeatures(
    src_file="van_der_pol_oscillator.cl", # This is your source file. 
    variable_names=["x", "y"],            # names for our variables
    parameter_names=["mu"],               # name for our parameters
    num_noise=1,                          # Number of Weiner noise variables
    observer=clode.Observer.threshold_2,  # Choose an observer
    stepper=clode.Stepper.dormand_prince, # Choose a stepper
    tspan=tspan,
)

# Define parameter values of interest (only a few for demonstration)
parameters = [-1, 0, 0.01, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0]

# array format as expected internally - see implementation details
pars_v = np.array([[par] for par in parameters])

# create initial conditions for each ODE instance
x0 = np.tile([1, 1], (len(parameters), 1))

# send the data to the OpenCL device
integrator.initialize(x0, pars_v)

# Run the simulation for tspan time, storing only the final state.
# Useful for integrating past transient behavior
integrator.transient()

# Continue the simulation, now measuring features of the solution
integrator.features()

# Get the results from the feature observer, print the period
observer_output = integrator.get_observer_results()
print(observer_output.get_var_mean("period"))
```

For more details, see the API reference [TODO]

## Implementation details

The Python library wraps a CPP library, clode_cpp_wrapper.[so|dll]
The CPP library assumes that the variables/parameters are grouped
by columns, i.e. if your variables are a, b and c,
the CPP library expects data in the format [aaaabbbbcccc].
The Python library expects data in the format
[[a, b, c], [a, b, c], [a, b, c], ...]

Note: There is a bug in the noise generation. It is currently not possible to
set the number of noise variables to 0. Set to 1 or more.
