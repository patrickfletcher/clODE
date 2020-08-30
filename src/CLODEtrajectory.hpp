/* clODE: a simulator class to run parallel ODE simulations on OpenCL capable hardware.
 * A clODE simulator solves on initial value problem over a grid of parameters and/or initial conditions. At each timestep,
 * "observer" rountine may be called to record/store/compute features of the solutions. Examples include storing the full
 * trajectory, recording the times and values of local extrema in a variable of the system, or directly computing other 
 * features of the trajectory.  
 */

//when compiling, be sure to provide the clODE root directory as a define:
// -DCLODE_ROOT="path/to/my/clODE/"

#ifndef CLODE_TRAJECTORY_HPP_
#define CLODE_TRAJECTORY_HPP_

#include "CLODE.hpp"
#include "clODE_struct_defs.cl"
#include "OpenCLResource.hpp"

#define __CL_ENABLE_EXCEPTIONS
#if defined(__APPLE__) || defined(__MACOSX)
#include "OpenCL/cl.hpp"
#else
#include <CL/cl.hpp>
#endif

#include <string>
#include <vector>

class CLODEtrajectory : public CLODE
{

protected:
    cl_int nStoreMax;
    std::vector<int> nStored;
    std::vector<double> t, x, dx, aux; //new result vectors
    size_t telements, xelements, auxelements;

    cl::Buffer d_t, d_x, d_dx, d_aux, d_nStored;
    cl::Kernel cl_trajectory;

    void initializeTrajectoryKernel();
    void resizeTrajectoryVariables(); //creates trajectory output global variables, called just before launching trajectory kernel

public:
    CLODEtrajectory(ProblemInfo prob, std::string stepper, bool clSinglePrecision, OpenCLResource opencl); //will construct the base class with same arguments
    ~CLODEtrajectory();

    //build program, set all problem data needed to run
    virtual void initialize(std::vector<double> newTspan, std::vector<double> newX0, std::vector<double> newPars, SolverParams<double> newSp);

    //simulation routines.
    //TODO: overload with newX0, newPars; all four?
    void trajectory(); //integrate forward an interval of duration (tf-t0)

    //Get functions
    std::vector<double> getT();
    std::vector<double> getX();
    std::vector<double> getDx();
    std::vector<double> getAux();
    std::vector<int> getNstored();
};

#endif //CLODE_HPP_
