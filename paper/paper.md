---
title: 'CLODE: A Python/C++ package for large-scale ODE simulation'
tags:
    - Python
    - C++
    - OpenCL
    - ODE
    - GPU
authors:
    - name: Patrick Fletcher
      affiliation: 1
    - name: Wolf Byttner
      affiliation: 2
affiliations:
    - name: National Institute of Health, USA
      index: 1
    - name: University of Exeter, UK
      index: 2
date: 2023-03-14
bibliography: paper.bib
---

# Summary

CLODE is a Python/C++ package for large-scale Ordinary
Differential Equation (ODE) simulation and analysis.
It uses OpenCL to parallelize ODE integration on a GPU, letting
users run over 100,000 simulations per second. The package is
designed to work with XPPAUT files, but can also be used as
a standalone solver with OpenCL-defined right-hand side equations.

# Statement of need

A large number of research questions concern ODEs, and scientists often
want to understand the aggregate influence of their parameters.
For meaningful statistical analysis, researchers need to run
millions to billions of simulations, which is not possible within
the time constraints of a typical research project.

Traditional ODE solvers in Python and Matlab, such as SciPy's
`odeint` and Matlab's `ode45`, are not designed for large-scale
simulations and take tens to hundreds of milliseconds to integrate
typical ODEs with up to a dozen state variables. CLODE is designed
to run on a GPU, and can solve ODEs in parallel, with each ODE
taking microseconds to integrate.

CLODE can also capture ODE features on the fly, without requiring
explicit trajectories, thus reducing memory requirements by over 99 %
for large-scale simulations with millions of ODEs. This lets researchers
run experiments on typical desktop hardware, without requiring dedicated
HPC clusters.

# Overview of the Library

CLODE contains two main components; an ODE feature observer, that
extracts features from the simulation on the fly, and a trajectory
simulator, that returns the explicit trajectory of a solution. The
feature observer is very memory efficient and can simulate [TODO FIND OUT HOW MANY] ODEs
simultaneously. The trajectory simulator is less memory efficient,
since it needs to retain the entire trajectory of each ODE, but it
is still much more performant than Scipy's `odeint`.

# Feature Observer

The feature observer can be initialised from Python and C++.
It supports several observers