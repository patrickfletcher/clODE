# clODE tests and benchmarks

## Tests verifying API function

There is a test file corresponding to each of the main python modules. These tests verify the function of the API.

## Tests for numerical correctness/accuracy

Verify that numerical results agree with target values to acceptable precision

### Tests for accuracy of solver and observer components

- Precision - single/mixed/double
- Steppers, adaptive step-size controllers
- Interpolation methods (for trajectory dense output, observers)
- Observer event detection methods (threshold crossings, local extrema)
- Observer online algorithms (running means, variances)

### Solutions of test problems

Test transient, trajectory, and features on test problems

- problems with analytic solutions
- gold-standard high-accuracy numerical solutions computed externally (i.e., with SciPy)

## Benchmarks measuring accuracy vs. performance

Running large ensemble simulations involves a trade-off between accuracy and computation time. The aim of these tests is to show how this trade-off depends on choices of methods and hyperparameters

- Follow SciML.ai: [Work-Precision Diagrams](https://docs.sciml.ai/SciMLBenchmarksOutput/stable/NonStiffODE/FitzhughNagumo_wpd/)
- Similarly for observers: online algorithms for measuring trajectory features

Other topics:

- floating point precision (single/compensated-single/double)
- time-stepping algorithms (fixed/adaptive, explicit/implicit, controllers)
- device performance, method hyperparamters for ensemble/time chunking (minibatching) depending on device/IVP
- memory layout/bandwidth (host-device, global/shared/registers), instruction throughput/streams, chunking/batching ensemble/timestepping, OpenCL NDRange, device dependent things
