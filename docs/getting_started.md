# Getting Started

## Basic concepts

clODE is a Python library for solving systems of ordinary differential
equations (ODEs) using OpenCL. It is designed to be easy to use, and
to provide a high-level interface for specifying ODEs and running simulations.

The primary use case for clODE is to simulate **ensembles** of initial value problems. This
is useful for studying the behavior of a system of ODEs as a function of
parameters, initial conditions, or other factors. By leveraging OpenCL, clODE enables significant speedups for this inherently parallel problem on any CPU, GPU, or other device with OpenCL support.

The clODE library includes three solvers:

- `Simulator`: This is the solver base class from which other solvers inherit. It supports advancing the ODE solution over a time interval, storing only the final state. This may be useful for studying long-term convergence behavior.
- `FeatureSimulator`: This solver measures features of the solution, such as the period of oscillation, on the fly, without storing the full trajectory data. The reduced storage can enable large ensembles to be rapidly evaluated.
- `TrajectorySimulator`: This solver stores the full trajectories of an ensemble of initial value problems.

clODE lets users specify ODEs using Python functions, and then converts these
functions to OpenCL code. This code is then executed on an OpenCL device, such
as a GPU or a multi-core CPU.

## Usage Example - computing the period of the Van der Pol oscillator

[The Van der Pol oscillator](https://en.wikipedia.org/wiki/Van_der_Pol_oscillator) can be written as a system of two differential equations:

$$
\dot{x} = y\\
\dot{y} = \mu(1-x^2)y - x
$$

Oscillations occur when $\mu>0$. Suppose we wish to measure the period of oscillations as $\mu$ varies.
First, we will need to implement the vector field above as function -
the right-hand-side (RHS) function - with the signature expected by clODE:

```python
def van_der_pol(float t,
                list[float] variables,
                list[float] parameters,
                list[float] derivatives,
                list[float] aux,
                list[float] weiner) -> None:


    # State variables
    x: float = variables[0]
    y: float = variables[1]

    # Parameters
    mu: float = parameters[0]

    # Differential equations
    dx: float = y
    dy: float = mu * (1 - x*x) * y - x

    # Differential outputs
    derivatives[0] = dx
    derivatives[1] = dy
```

Note that the Python function must be **fully typed** and
must have this exact signature. This function will be converted to OpenCL and written
to a file called `clode_rhs.cl`. This file will then be loaded
into an OpenCL program and executed on the OpenCL device.

Note that this differs from the signature expected by Scipy's ```solve_ivp```;
[TODO] we provide a wrapper to support using the same vector field function in both Scipy and clODE.

In the above case, the output OpenCL function would look like this:

```c
void getRHS(const realtype t,
            const realtype variables[],
            const realtype parameters[],
            realtype derivatives[],
            realtype aux[],
            const realtype wiener[]) {

    /* State variables */
    realtype x = variables[0];
    realtype y = variables[1];

    /* Parameters */
    realtype mu = parameters[0];

    /* Differential equations */
    realtype dx = y;
    realtype dy = mu * (1 - x*x) * y - x;

    /* Differential outputs */
    derivatives[0] = dx;
    derivatives[1] = dy;
}
```

Note that this is a simple C-language function.
The ```realtype``` type declaration is a macro that will expand to
```float``` or ```double```, depending on whether clODE is configured for
single or double precision floats. clODE defaults to single precision.

Note - clODE supports many additional use cases. For reference, see
[Specifying systems of ODEs](specifying_odes.md).

[//]: # (This function supports additional use cases not used in this example: ```t``` for time-dependent ODE terms &#40;non-autonomous systems&#41;, ```aux``` for auxiliary readout variables, and stochastic terms via ```w```, which provides Wiener variables &#40;$w \sim Normal \; &#40;0,dt&#41;$&#41;. If needed, one can also declare additional plain C-languange functions in the same file, preceding ```getRHS``, and use them inside getRHS.)

Next we will use a python script to define our parameters and set up the numerical simulation. Here we use clODE's feature detection mode - several features of the ODE solution, including the period of oscillation, will be measured "on the fly", without storing the trajectory itself.

```python
from clode import Stepper, Observer, CLODEFeatures
import numpy as np

# time span for our simulation
tspan = (0.0, 1000.0)

# Create the clODE feature extractor
integrator = CLODEFeatures(
    src_file="van_der_pol_oscillator.cl",  # This is your source file. 
    variable_names=["x", "y"],  # names for our variables
    parameter_names=["mu"],  # name for our parameters
    observer=Observer.threshold_2,  # Choose an observer
    stepper=Stepper.rk4,  # Choose a stepper
    tspan=tspan,
)

# Define parameter values of interest (only a few for demonstration)
mu = [0.01, 0.5, 2.0, 4.0]

# array format as expected internally - see implementation details
P0 = np.array([[u] for u in mu])

# create initial conditions for each ODE instance
x0 = np.tile([1, 1], (len(mu), 1))

# send the data to the OpenCL device
integrator.set_ensemble(x0, P0)

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

## Trajectories

We can also compute and store full trajectories in parallel. This requires significantly more memory, though, but is important for validating the above feature results.  Continuing the example above, we next compute and plot the trajectories for the four parameters specified.

```python
from clode import CLODETrajectory
import matplotlib.plt as plt

# Create the clODE trajectory solver
integrator = CLODETrajectory(
    src_file="van_der_pol_oscillator.cl",  # This is your source file. 
    variable_names=["x", "y"],  # names for our variables
    parameter_names=["mu"],  # name for our parameters
    stepper=Stepper.rk4,  # Choose a stepper
    tspan=tspan,
)

# send the data to the OpenCL device
integrator.set_ensemble(x0, P0)

# Run the simulation for tspan time, storing only the final state.
# Useful for integrating past transient behavior
integrator.transient()

# Continue the simulation, now storing the trajectories
trajectories = integrator.trajectory()

# plot
fig, ax = plt.subplots(4, 1, sharex=True, sharey=True)

for i, trajectory in enumerate(trajectories):
    ax[i].plot(trajectory["t"], trajectory.x['x'])

ax[1].set_ylabel('x')
ax[-1].set_xlabel('time')
plt.show()
```

## Implementation details

The Python library wraps a CPP library, clode_cpp_wrapper.[so|dll]
The CPP library assumes that the variables/parameters are grouped
by columns, i.e. if your variables are a, b and c,
the CPP library expects data in the format [aaaabbbbcccc].
The Python library expects data in the format
[[a, b, c], [a, b, c], [a, b, c], ...]
