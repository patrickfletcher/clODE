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

clODE is an efficient computational tool designed for parallel solving of ordinary differential equation (ODE) ensembles using OpenCL. It lets users define their ODE system and the ensemble of parameter sets and initial conditions in Python. By leveraging OpenCL, significant speedups can be obtained for this inherently parallel problem on CPUs, GPUs, and other OpenCL-capable devices.

The public Python API is built around three simulator classes:

- `Simulator` advances an ensemble and keeps only the final state.
- `FeatureSimulator` computes trajectory features on the fly through a stateful observer, without storing the full trajectory.
- `TrajectorySimulator` stores full trajectory samples.

clODE offers flexibility in simulator deployment across different hardware, allowing, for example, the `FeatureSimulator` to operate on a GPU while the `TrajectorySimulator` runs on a CPU.

The repository is currently transitioning from a legacy C++/Bazel host runtime to a Python-owned PyOpenCL backend. The public API remains stable during that migration. Today, the default backend is still the legacy path, while the optional `clode[pyopencl]` dependency enables the transition backend for contributor workflows and validation.

## Installation

See [installation](https://patrickfletcher.github.io/clODE/install/) for current installation instructions.

## Getting Started

See [Getting Started](https://patrickfletcher.github.io/clODE/getting_started/) for current Python examples.

## Source

The source code is available on [GitHub](https://github.com/patrickfletcher/clODE).
