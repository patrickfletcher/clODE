# clODE - an OpenCL based tool for solving ordinary differential equations (ODEs)

[![Python](https://img.shields.io/pypi/pyversions/clode.svg)](https://badge.fury.io/py/clode)
[![PyPI version](https://badge.fury.io/py/clode.svg)](https://badge.fury.io/py/clode)
[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/patrickfletcher/clODE/badge)](https://securityscorecards.dev/viewer/?uri=github.com/patrickfletcher/clODE)
![Windows](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_build_windows.yml/badge.svg)
![Mac](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_test_mac.yml/badge.svg)
![Linux](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_build_linux.yml/badge.svg)

**`Documentation`** |
------------------- |
[![Documentation](https://img.shields.io/badge/api-reference-blue.svg)](https://patrickfletcher.github.io/clODE/) |


clODE is an efficient computational tool designed for parallel solving
of ordinary differential equation (ODE) ensembles using OpenCL.
It lets users define their ODE system and the ensemble of parameter sets and initial conditions in Python.  By leveraging OpenCL, significant speedups can be obtained for this inherently parallel problem on any CPU, GPU, or other device with OpenCL support. Two primary modes of simulation are supported:

- FeatureSimulator computes features of ODE trajectories, such as oscillation period, on-the-fly, without storing the trajectory data, facilitating extensive parameter analyses with considerable computational speed improvements.
- TrajectorySimulator stores the full trajectory data.

clODE offers flexibility in simulator deployment across different hardware,
allowing, for example, the FeatureSimulator to operate on a GPU while the
TrajectorySimulator runs on a CPU.

Developed in C++ and OpenCL, clODE is accessible for direct use in C++
applications or through a Python interface. The library compiles with bazel
and bazelisk, and works on Linux, Windows, and MacOS platforms.

## Installation

See [installation](https://patrickfletcher.github.io/clODE/install/) for instructions on how to install CLODE.

## Getting Started

See [Getting Started](https://patrickfletcher.github.io/clODE/getting_started/) for an example of clODE usage.

## Source

The source code is available on [GitHub](https://github.com/patrickfletcher/clODE).

