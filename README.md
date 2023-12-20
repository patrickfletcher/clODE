# clODE - an OpenCL based tool for solving ordinary differential equations (ODEs)

[![Python](https://img.shields.io/pypi/pyversions/clode.svg)](https://badge.fury.io/py/clode)
[![PyPI](https://badge.fury.io/py/clode.svg)](https://badge.fury.io/py/tensorflow)
[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/patrickfletcher/clODE/badge)](https://securityscorecards.dev/viewer/?uri=github.com/patrickfletcher/clODE)
![Windows](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_build_windows.yml/badge.svg)
![Mac](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_test_mac.yml/badge.svg) 
![Linux](https://github.com/patrickfletcher/clODE/actions/workflows/bazel_build_linux.yml/badge.svg)

clODE is a tool for solving ordinary differential equations (ODEs) using OpenCL. It is tailored to numerically solving many instances of a given ODE system in parallel, each with different parameter sets and/or initial conditions.

The ODE solver runs entirely on the OpenCL device, supporting independent solver state per simulation (e.g., adaptive timesteps). clODE can return the full trajectories, though this is somewhat memory intensive. Alternatively, clODE supports computing features of the ODE trajectory (e.g., oscillation period) on the fly without storing the trajectory itself, enabling much larger parameter sweeps to be run with significant speedup over serial computation.

clODE is written in C++ and OpenCL, and can be used directly in C++ programs or via the provided Python interface. The library is compiled using bazel and bazelisk, and it runs on Linux, Windows and MacOS.

## Installation

See [installation](docs/install.md) for instructions on how to install CLODE.

## Getting Started

See [Getting Started](docs/getting_started.md) for an example of clODE usage.

## Source

The source code is available on [GitHub](https://github.com/patrickfletcher/clODE).
