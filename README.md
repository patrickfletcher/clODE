# clODE

[![Python](https://img.shields.io/pypi/pyversions/clode.svg)](https://badge.fury.io/py/clode)
[![PyPI version](https://badge.fury.io/py/clode.svg)](https://badge.fury.io/py/clode)
[![CI](https://github.com/patrickfletcher/clODE/actions/workflows/ci.yml/badge.svg)](https://github.com/patrickfletcher/clODE/actions/workflows/ci.yml)
[![Docs](https://github.com/patrickfletcher/clODE/actions/workflows/deploy_mkdocs_pages.yml/badge.svg)](https://github.com/patrickfletcher/clODE/actions/workflows/deploy_mkdocs_pages.yml)
[![Release](https://github.com/patrickfletcher/clODE/actions/workflows/python_package_release.yml/badge.svg)](https://github.com/patrickfletcher/clODE/actions/workflows/python_package_release.yml)
[![OpenSSF Scorecard](https://api.securityscorecards.dev/projects/github.com/patrickfletcher/clODE/badge)](https://securityscorecards.dev/viewer/?uri=github.com/patrickfletcher/clODE)

clODE is a Python package for large-scale simulation of ordinary differential equation ensembles on OpenCL-capable CPUs and GPUs. It is built for workloads where you want to sweep parameters, run many independent systems in parallel, and choose between final-state simulation, online feature extraction, or stored trajectories.

clODE supports:

- typed Python RHS functions, OpenCL source files, and XPP models
- deterministic and stochastic ODE systems
- feature extraction without storing full trajectories
- explicit OpenCL device inspection and runtime selection from Python

## Installation

```bash
pip install clode
```

An OpenCL runtime for your target device is required. See the [installation guide](https://patrickfletcher.github.io/clODE/install/) for platform notes, verification steps, and source-install details.

## Quick Start

The example below measures oscillation period for an ensemble of Van der Pol systems without storing full trajectories.

```python
from typing import List

import clode
import numpy as np


def van_der_pol(
    t: float,
    variables: List[float],
    parameters: List[float],
    derivatives: List[float],
    aux: List[float],
    wiener: List[float],
) -> None:
    x: float = variables[0]
    y: float = variables[1]
    mu: float = parameters[0]

    derivatives[0] = y
    derivatives[1] = mu * (1.0 - x * x) * y - x


integrator = clode.FeatureSimulator(
    rhs_equation=van_der_pol,
    variables={"x": 1.0, "y": 1.0},
    parameters={"mu": 0.1},
    observer=clode.Observer.threshold_2,
    stepper=clode.Stepper.dormand_prince,
    t_span=(0.0, 1000.0),
)

integrator.set_ensemble(parameters={"mu": np.array([0.01, 0.5, 2.0, 4.0])})
integrator.transient()

observer_output = integrator.features()
print(observer_output.get_var_mean("period"))
```

For a fuller walkthrough, see [Getting started](https://patrickfletcher.github.io/clODE/getting_started/).

## Documentation

- [Documentation site](https://patrickfletcher.github.io/clODE/)
- [Installation guide](https://patrickfletcher.github.io/clODE/install/)
- [Getting started](https://patrickfletcher.github.io/clODE/getting_started/)
- [Examples](https://patrickfletcher.github.io/clODE/examples/)
- [Performance notes](https://patrickfletcher.github.io/clODE/performance_notes/)
- [Feature extraction](https://patrickfletcher.github.io/clODE/feature_extraction/)
- [Trajectory simulation](https://patrickfletcher.github.io/clODE/trajectory_simulation/)
- [Specifying ODE systems](https://patrickfletcher.github.io/clODE/specifying_odes/)
- [API reference](https://patrickfletcher.github.io/clODE/api_reference/)

## Repository Layout

- [clode/](https://github.com/patrickfletcher/clODE/tree/main/clode) contains the maintained Python package source.
- [docs/](https://github.com/patrickfletcher/clODE/tree/main/docs) contains the MkDocs documentation source.
- [examples/](https://github.com/patrickfletcher/clODE/tree/main/examples) contains runnable examples and sample model files.
- [test/](https://github.com/patrickfletcher/clODE/tree/main/test) contains the regression suite.
- [paper/](https://github.com/patrickfletcher/clODE/tree/main/paper) contains the software paper materials.

## Contributing

See [CONTRIBUTING.md](https://github.com/patrickfletcher/clODE/blob/main/CONTRIBUTING.md) for local setup, test bundles, docs commands, and contributor expectations.

## License

clODE is distributed under the MIT License. See [LICENSE](https://patrickfletcher.github.io/clODE/LICENSE/).
