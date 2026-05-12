---
title: 'clODE: A Python package for large-scale ensemble ODE simulation and online feature extraction with OpenCL'
tags:
  - Python
  - OpenCL
  - ODE
  - GPU
authors:
  - name: Patrick Fletcher
    affiliation: 1
  - name: Wolf Byttner
    affiliation: "2, 3"
affiliations:
  - name: National Institute of Health, United States of America
    index: 1
  - name: University of Exeter, United Kingdom
    index: 2
  - name: Planar AI, United Kingdom
    index: 3
date: 2023-03-14
bibliography: paper.bib
---

<!-- markdownlint-disable MD025 -->

# Summary

clODE is a Python package for simulating large ensembles of ordinary
differential equation (ODE) systems on OpenCL-capable CPUs and GPUs.
It is built around three related workflows: advancing an ensemble while
retaining only the final state, computing features and events on the device
during integration, and storing full trajectories for later inspection.
Models can be supplied as typed Python right-hand-side functions, OpenCL
source files, or XPP models.

The package targets workloads where many closely related solves must be run
across parameter grids or repeated ensembles, especially when users care more
about aggregate behavior than about retaining every time sample. clODE was
developed to study populations of ODE models in computational neuroscience
[@fletcher:2016; @fletcher:2017] and remains focused on high-throughput,
feature-oriented ensemble simulation from Python.

# Statement of need

Many scientific modelling problems are not limited by the accuracy of a single
ODE solve; they are limited by the need to run large ensembles of related solves
across parameter sweeps, initial-condition ensembles, or repeated stochastic
realizations. In these settings, storing every trajectory sample is often
unnecessary and can dominate both memory use and runtime.

clODE addresses this workflow directly. It exposes an OpenCL-backed ensemble
simulation engine through a Python package and supports on-device feature and
event extraction so users can collect periods, extrema, counts, and event data
without always storing full trajectories. The package is designed to keep model
authoring approachable by accepting typed Python right-hand-side functions,
OpenCL source, and XPP models while still targeting OpenCL-capable accelerators
and CPUs.

This combination is useful for researchers who need a workflow-oriented tool for
large ensemble simulation rather than a maximally broad general-purpose ODE
ecosystem. The package is especially aimed at studies where model structure is
fixed, the number of solves is large, and the desired outputs are summary
statistics or event-based features rather than dense trajectories.

## Alternatives considered

SciPy's `odeint` and `solve_ivp` functions [@2020SciPy-NMeth] are strong
general-purpose tools for individual solves and small to medium parameter scans,
but they do not target OpenCL-backed ensemble execution or on-device feature
extraction. DifferentialEquations.jl and Python access paths such as diffeqpy
[@DifferentialEquations.jl-2017] offer a much broader solver ecosystem, but with
a different deployment and ergonomics story than a focused Python package built
around OpenCL ensemble workflows.

Machine-learning-oriented ODE packages such as torchdiffeq
[@chen2018neuralode] and torchode [@lienen2022torchode] emphasize differentiable
or tensor-native batched solves inside deep-learning stacks. That is a different
use case from clODE's observer- and feature-oriented ensemble workflows. XPPAUT,
meanwhile, remains an important authoring and analysis environment for dynamical
systems models, but it is not itself an OpenCL-backed engine for large ensemble
simulation.

clODE therefore occupies a narrower niche: a Python package for high-throughput
ensemble simulation on OpenCL devices, with support for model ingestion from
Python, OpenCL source, and XPP, and with built-in attention to feature extraction
without mandatory full-trajectory storage.

# Current applications and future work

clODE is currently being used in studies of pituitary-cell dynamics and related
computational neuroscience workflows where large ensembles of closely related ODE
models must be explored [@fletcher:2016; @fletcher:2017]. The current package
also supports XPP-based workflows, making it possible to connect existing model
authoring practices to OpenCL-backed ensemble execution.

Near-term work focuses on strengthening reproducibility and public evaluation:
clearer benchmark material, tighter documentation of workflow tradeoffs, and
continued refinement of observer- and continuation-related semantics so the
public story remains aligned with the supported implementation.

# Acknowledgements

We would like to thank Joel Tabak and Richard Bertram for extensively testing
CLODE. We would also like to thank the reviewers for their helpful comments.

# References
