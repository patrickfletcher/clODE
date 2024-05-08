//
// Created by Patrick Fletcher 2017
//

#ifndef CLODE_TRAJECTORY_HPP_
#define CLODE_TRAJECTORY_HPP_

#include "CLODE.hpp"
#include "clODE_struct_defs.cl"
#include "OpenCLResource.hpp"

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include "OpenCL/cl2.hpp"

#include <string>
#include <vector>

class CLODEtrajectory : public CLODE
{

protected:
    cl_int nStoreMax;
    std::vector<cl_int> nStored;
    std::vector<cl_double> t, x, dx, aux; // new result vectors
    size_t telements, xelements, auxelements;

    cl::Buffer d_t, d_x, d_dx, d_aux, d_nStored;
    cl::Kernel cl_trajectory;

    void resizeTrajectoryVariables(); // creates trajectory output global variables, called just before launching trajectory kernel

public:
    CLODEtrajectory(ProblemInfo prob, std::string stepper, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot); // will construct the base class with same arguments
    CLODEtrajectory(ProblemInfo prob, std::string stepper, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot);
    virtual ~CLODEtrajectory();

    void buildCL() override; // build program and create kernel objects

    // simulation routine. TODO: overloads?
    void trajectory(); // integrate forward an interval of duration (tf-t0)

    // Get functions
    std::vector<cl_double> getT();
    std::vector<cl_double> getX();
    std::vector<cl_double> getDx();
    std::vector<cl_double> getAux();
    std::vector<cl_int> getNstored();
};

#endif // CLODE_HPP_
