# clODE - an OpenCL based tool for solving ordinary differential equations (ODEs)

clODE is a tool for solving ordinary differential equations (ODEs) using OpenCL.
It is tailored to numerically solving a given ODE system for large collections of parameter sets and/or initial conditions in parallel. The ODE solver runs entirely on the OpenCL device, supporting independent solver state per simulation (e.g., adaptive timesteps). clODE also supports computing features of the ODE trajectory (e.g., oscillation period) on the fly without storing the trajectory itself, enabling very large parameter sweeps to be run with significant speedup over serial computation.

clODE is written in C++ and OpenCL, and has a Python interface.
The ODE system's vector field function is written in OpenCL,
and supports exporting auxiliary/readout variables. The library is compiled
using bazel and bazelisk, and it runs on Linux, Windows and MacOS.

## Source

The source code is available on [GitHub](https://github.com/patrickfletcher/clODE).
