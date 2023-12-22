# clODE Matlab interface

## Note: currently unmaintained

The Matlab interface will migrate to direct calls to the python version from Matlab.

## Overview

clODE solves many instances of the same ODE system given different input paramter values and/or initial conditions.

This document describes the Matlab MEX-file interface provided for running clODE from within Matlab.

## Installation

Download/clone this repository to a folder of your choice. Add the `/clODE/matlab` folder to your matlab path.

To compile the mex files, first edit the appropriate path variables in `compileCLODEmex.m` to point to your OpenCL installation:

``` matlab
  opencl_include_dir = '/path/to/cl.hpp'; %all platforms
  opencl_lib_dir = '/path/to/libOpenCL.so'; %Linux
  opencl_lib_dir = '/path/to/OpenCL.lib'; %Windows
```

Next, run the script. Requires a [C++ compiler compatible with Matlab MEX-file compilation](https://www.mathworks.com/support/requirements/supported-compilers.html).

The specification of ODE systems can be done using XPPAUT style ODE files. To use this feature, please also [download XPPToolbox](https://github.com/patrickfletcher/xppToolbox) and [add it to your Matlab path](https://www.mathworks.com/help/matlab/ref/addpath.html).

## Basics

ODEs can be integrated with two main modes of result storage. Trajectories can be recorded in full (using the `clODEtrajectory` class), which may result in large memory requirements if the number of IVP instances gets large. Alternatively, features of the trajectories can be computed on the fly during integration without storing the full trajectories (using the `clODEfeatures` class).

These classes can be used in scripts or via graphical interface elements, as described in the following sections. Some simple test examples are provided in the `/clODE/samples/` subdirectory.

## Setting up a solver

clODE objects are created via the class constructors. At minimum, the constructors require one input, pointing to the ODE system definition. Additionally, the numerical precision, OpenCL device, ODE time stepping algorithm, and feature-detection method (for clODEfeatures only) can be specified. Examples:

``` matlab
clo=clODE(odefile); %using an XPP ODE-file to specify the ODE, using default solver specification: precision='single', selectedDevice=1, stepper='dorpri5'.

openclDevices=queryOpenCL(); %array of structs describing OpenCL devices available
precision='double'; %'single' {default} or 'double' 
selectedDevice=2; %1-based index into openclDevices.
stepper='rk4'; %a valid time stepping algorithm name 
clo=clODE(odefile,precision,selectedDevice,stepper); 

clo=clODE(odestruct); %using a struct to specify the ODE system, e.g. output of ode2cl
```

The returned `clo` object is used to run simulations. ODE system information is stored in the struct `clo.prob`.

### ODE system specification

The easiest way to specify an ODE system is to use an ODE file as needed for XPPAUT. The `ode2cl` function from [XPPToolbox](https://github.com/patrickfletcher/xppToolbox) to parse XPPAUT files and automatically generate the OpenCL code required for execution on OpenCL devices. It also returns a struct containing detailed information about the ODE system for use in your programs that can be used to initialize a clODE solver.

Handwritten OpenCL files are also supported. At minimum a struct containing the following information is used to specify the ODE system:

``` matlab
problem.clRHSfilename='/path/to/RHSfile.cl'; %An OpenCL function that computes the system's vector field
problem.nVar=integer; %number of state variables
problem.nPar=integer; %number of parameters that may vary
problem.nAux=integer; %number of auxiliary output quantities (user defined functions of state variables)
problem.nWiener=integer; %number of Wiener random variables for stochastic simulations
```

The `RHSfile.cl` file must implement a function with the following signature:

``` c
void getRHS(realtype t, realtype x_[], realtype p_[], realtype dx_[], realtype aux_[], realtype w_[]);
// realtype is float or double
// Input arguments: x_[] - state variable array, p_[] - parameter array, w_[] - array of Wiener variable values provided by clODE if nWiener>0. 
// Output arguments: dx_[] - slopes, aux_[] - values of auxiliary output quantities
```

## Setting up problem data

## clODEtrajectory

## clODEfeatures

The following

###
