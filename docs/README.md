# clODE - an OpenCL based tool for solving ordinary differential equations (ODEs)

[![Python](https://img.shields.io/pypi/pyversions/clode.svg)](https://badge.fury.io/py/clode)
[![PyPI version](https://badge.fury.io/py/clode.svg)](https://badge.fury.io/py/clode)
[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/patrickfletcher/clODE/badge)](https://securityscorecards.dev/viewer/?uri=github.com/patrickfletcher/clODE)
![Windows](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_build_windows.yml/badge.svg)
![Mac](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_test_mac.yml/badge.svg)
![Linux](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_build_linux.yml/badge.svg)

clODE is an efficient computational tool designed for parallel solving
of ordinary differential equation (ODE) ensembles by leveraging OpenCL.
It lets users define ODE systems in Python and enhances numerical
simulation speeds up to 500 times faster than scipy's solve_ivp by
utilising graphics cards, on-the-fly feature extraction
and parallel processing.

The software features two primary simulators: the FeatureSimulator,
which analyses ODE trajectory characteristics (like oscillation periods)
in real-time without needing to store the trajectory, facilitating extensive
parameter analyses with considerable computational speed improvements.
Conversely, the TrajectorySimulator provides detailed trajectory data.
clODE offers flexibility in simulator deployment across different hardware,
allowing, for example, the FeatureSimulator to operate on a GPU while the
TrajectorySimulator runs on a CPU.

Developed in C++ and OpenCL, clODE is accessible for direct use in C++
applications or through a Python interface. The library compiles with bazel
and bazelisk, and works on Linux, Windows, and MacOS platforms.

## Installation

See [installation](install.md) for instructions on how to install CLODE.

## Getting Started

See [Getting Started](getting_started.md) for an example of clODE usage.

## Source

The source code is available on [GitHub](https://github.com/patrickfletcher/clODE).

