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

clODE is a Python package for large-scale ensemble simulation of ordinary differential equations (ODEs) on OpenCL-capable CPUs and GPUs. It is built for parameter sweeps and ensemble studies where you want to run many closely related ODE solves in parallel and extract summary features or event data on the device, avoiding the memory and runtime cost of storing every trajectory sample.

The package works through three complementary workflows: (1) advancing an ensemble and retaining only final states for convergence studies, (2) computing features and events on the device during integration through stateful observers that summarize trajectories in situ, and (3) storing full trajectories for visualization and post-analysis when trajectory data is needed.

Models can be supplied as typed Python right-hand-side functions, OpenCL source files, or XPP models. The observer architecture includes seven canonical families—summary statistics, threshold crossing (absolute and warmup-normalized), Schmitt triggering, local extrema, and neighborhood-return triggers—designed to cover common dynamical-systems workflows without storing full trajectories.

The package also emphasizes numerically robust single-precision ensemble workflows by keeping solver time and on-device summaries in compensated relative forms where that materially improves long OpenCL runs.

# Statement of Need

Large-scale ensemble simulation of dynamical systems is a common research need: parameter sweeps, bifurcation analysis, sensitivity studies, and repeated stochastic realizations all require running many related ODE solves. However, storing trajectories for every ensemble member is often unnecessary and can dominate both memory and runtime, especially when the research question is about ensemble behavior or specific features such as oscillation periods, local extrema, or event detection rather than full trajectory inspection.

clODE addresses this directly through **on-device feature extraction**. Instead of computing trajectories and storing them to the host, clODE runs stateful observer kernels on the GPU or accelerator that accumulate features during integration—collecting periods, extrema, event counts, and event timestamps without keeping trajectory samples in device memory. This combination reduces memory footprint, avoids GPU-to-host I/O bottlenecks, and lets users scale to much larger ensembles than trajectory-focused approaches allow.

The package keeps model authoring accessible by accepting typed Python RHS functions, OpenCL source files, and XPP models, while maintaining a PyOpenCL-only runtime for broad hardware support without vendor lock-in. It exposes an explicit observer architecture through a semantic config layer so users can choose the right observer family for their workflow—threshold crossing, Schmitt triggering, local extrema detection, or full state summaries—without rewriting simulation code.

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

# Current Applications and Future Work

clODE is currently in use for computational neuroscience workflows including pituitary-cell population studies [@fletcher:2016; @fletcher:2017], bifurcation parameter sweeps, and phase-response analysis where large ensembles of closely related ODE models must be explored efficiently. The package's support for XPP model ingestion makes it a natural bridge between XPPAUT-based model authoring and OpenCL-backed ensemble execution on modern hardware.

Near-term development focuses on strengthening the public evidence base and user experience. Key priorities are: benchmark and reproducibility material that demonstrates speedups on real research workloads; clearer numerical-accuracy and observer-tradeoff material around the compensated single-precision paths already used in the package; and expanded guidance for common dynamical-systems workflows such as bifurcation analysis and phase-space exploration.

The current architecture deliberately stays narrow—high-throughput ensemble simulation with on-device features—allowing it to stay maintainable and numerically robust rather than attempting to be a universal ODE solver platform.

# Acknowledgements

We would like to thank Joel Tabak and Richard Bertram for extensively testing
CLODE. We would also like to thank the reviewers for their helpful comments.

# References
