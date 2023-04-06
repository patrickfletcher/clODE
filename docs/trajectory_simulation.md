# Trajectory simulation

CLODE can simulate ODE trajectories using the CLODETrajectory class.

## Example

The following example simulates the Van der Pol oscillator using the dormand_prince45 integrator.

### OpenCL RHS function

```c
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

### Python

```python
import clode
import numpy as np

def ornl_thompson_a1():
    num_simulations = 1

    m = 1 / 4
    w = 8
    k = 2
    H = 10

    tspan = (0.0, H / 2.)

    integrator = clode.CLODETrajectory(
        src_file="test/ornl_thompson_a1.cl",
        variable_names=["y1", "y2"],
        parameter_names=["m", "w", "k", "H"],
        aux=["g1"],
        num_noise=0,
        dt=0.001,
        dtmax=0.001,
        stepper=clode.Stepper.rk4,
        tspan=tspan,
    )

    parameters = [m, w, k, H]

    x0 = np.tile([0, 0], (num_simulations, 1))
    pars_v = np.tile(parameters, (num_simulations, 1))

    integrator.initialize(x0, pars_v)

    integrator.trajectory()

    time_steps = integrator.get_time_steps()
    return integrator.get_trajectory()
```
