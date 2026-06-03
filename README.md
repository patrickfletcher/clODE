# clODE

[![Python](https://img.shields.io/pypi/pyversions/clode.svg)](https://badge.fury.io/py/clode)
[![PyPI version](https://badge.fury.io/py/clode.svg)](https://badge.fury.io/py/clode)
[![CI](https://github.com/patrickfletcher/clODE/actions/workflows/ci.yml/badge.svg)](https://github.com/patrickfletcher/clODE/actions/workflows/ci.yml)
[![Docs](https://github.com/patrickfletcher/clODE/actions/workflows/deploy_mkdocs_pages.yml/badge.svg)](https://github.com/patrickfletcher/clODE/actions/workflows/deploy_mkdocs_pages.yml)
[![Release](https://github.com/patrickfletcher/clODE/actions/workflows/python_package_release.yml/badge.svg)](https://github.com/patrickfletcher/clODE/actions/workflows/python_package_release.yml)
[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/patrickfletcher/clODE/badge)](https://securityscorecards.dev/viewer/?uri=github.com/patrickfletcher/clODE)

clODE is a Python package for large-scale simulation of ODE ensembles on OpenCL-capable CPUs and GPUs. It is built for parameter sweeps and ensemble studies where you need to run many related solves in parallel and extract features or summary statistics on the device, avoiding the memory and I/O cost of storing full trajectories. The package also includes numerically hardened single-precision timekeeping and summary accumulation paths for long ensemble workflows, with the tradeoffs documented in the numerical-accuracy guide.

clODE is particularly useful when you want to:

- sweep parameters across a dense grid and study ensemble behavior rather than individual trajectories
- compute periods, extrema, event counts, or event timestamps without storing all time samples
- leverage OpenCL for high throughput while keeping model authoring in Python

The package supports:

- models defined as typed Python RHS functions, OpenCL source files, or XPP files
- deterministic and stochastic ODE systems
- on-device feature and event extraction through stateful observers
- explicit OpenCL device inspection and single-device runtime selection from Python

## Installation

```bash
pip install clode
```

An OpenCL runtime for your target device is required. See the [installation guide](https://patrickfletcher.github.io/clODE/install/) for platform notes, verification steps, and source-install details.

## Quick Start

Here is a minimal example: create an ensemble of Van der Pol oscillators with different damping parameters, run an on-device summary pass without storing trajectories, and report the window-mean `x` value for each ensemble member:

```python
from typing import List
import clode
import numpy as np

def van_der_pol(t, variables, parameters, derivatives, aux, wiener):
    x, y = variables[0], variables[1]
    mu = parameters[0]
    derivatives[0] = y
    derivatives[1] = mu * (1.0 - x*x) * y - x

simulator = clode.FeatureSimulator(
    rhs_equation=van_der_pol,
    variables={"x": 1.0, "y": 1.0},
    parameters={"mu": 0.1},
    t_span=(0.0, 200.0),
)

simulator.set_ensemble(parameters={"mu": np.array([0.01, 0.5, 2.0, 4.0])})
summary = simulator.features()
print(summary.get_var_mean("x"))
```

For a step-by-step walkthrough, richer event observers, and multiple workflow examples, see [Getting started](https://patrickfletcher.github.io/clODE/getting_started/) and [Feature extraction](https://patrickfletcher.github.io/clODE/feature_extraction/).

## Documentation

Start with [Installation](https://patrickfletcher.github.io/clODE/install/) to set up your OpenCL runtime, then [Getting started](https://patrickfletcher.github.io/clODE/getting_started/) for your first end-to-end example.

Workflow guides:

- [Feature extraction](https://patrickfletcher.github.io/clODE/feature_extraction/) — compute periods, extrema, and event data without storing trajectories
- [Trajectory simulation](https://patrickfletcher.github.io/clODE/trajectory_simulation/) — store time samples for plotting and inspection
- [Continuation and repeated runs](https://patrickfletcher.github.io/clODE/continuation/) — default continued solves, explicit resets, and split-window behavior
- [Specifying ODE systems](https://patrickfletcher.github.io/clODE/specifying_odes/) — Python, OpenCL, and XPP model definition
- [Stochastic simulation](https://patrickfletcher.github.io/clODE/Ornstein-Uhlenbeck_process/) — add Wiener-process terms

Reference and examples:

- [API reference](https://patrickfletcher.github.io/clODE/api_reference/)
- [Examples](https://patrickfletcher.github.io/clODE/examples/)
- [Numerical accuracy](https://patrickfletcher.github.io/clODE/numerical_accuracy/)
- [Performance notes](https://patrickfletcher.github.io/clODE/performance_notes/)
- [Full docs](https://patrickfletcher.github.io/clODE/)

## Repository Layout

- `clode/` — maintained Python package source.
- `docs/` — MkDocs documentation.
- `examples/` — runnable workflow examples.
- `test/` — regression suite, including `core_numerics/` for the release gate.
- `paper/` — software paper source.

## Contributing

**Contributing:**

See [CONTRIBUTING.md](https://github.com/patrickfletcher/clODE/blob/main/CONTRIBUTING.md) for local setup, test bundles, docs commands, and contributor expectations. The package uses a release gate over `test/core_numerics/`, an extended test suite, and automated CI for reproducibility.

## License

clODE is distributed under the MIT License. See [LICENSE](https://patrickfletcher.github.io/clODE/LICENSE/).
