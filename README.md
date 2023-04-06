# clODE - an OpenCL based tool for solving ordinary differential equations (ODEs)

clODE is a tool for solving ordinary differential equations (ODEs) using OpenCL.
It lets users simulate 10,000s of ODEs simultaneously on a GPU, 
and generates data three orders of magnitude faster than Matlab's ODE solvers
or scipy's odeint.

clODE is written in C++ and OpenCL, and has a Python interface.
The right-hand-side (RHS) function is written in OpenCL,
and is relatively simple in structure. The library is compiled
using bazel and bazelisk, and it runs on Linux, Windows and MacOS.

## Source

The source code is available on [GitHub](https://github.com/patrickfletcher/clODE).