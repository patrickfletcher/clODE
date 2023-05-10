/* clODE: a simulator class to run parallel ODE simulations on OpenCL capable hardware.
 * A clODE simulator solves an initial value problem for a set of (parameters, initial conditions). At each timestep,
 * "observer" rountine may be called to record/store/compute features of the solutions. Examples include storing the full
 * trajectory, recording the times and values of local extrema in a variable of the system, or directly computing other 
 * features of the trajectory.  
 */

//TODO: namespaces?

//TODO: break up computation (timespan) into chunks that don't crash the system. some fine-grained max time chunk to run a kernel, a while loop in C++ to do all chunks

//TODO: choosing specific RNG - get nRNGstate using a switch

//TODO: device-to-device transfers instead of overwriting x0?

//TODO: separate flags for initialized state and built state

//TODO: use kernel.getWorkGroupInfo with CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE to offer the user a hint of nPts to use

//when compiling, be sure to provide the clODE root directory as a define:
// -DCLODE_ROOT="path/to/my/clODE/"

#ifndef CLODE_HPP_
#define CLODE_HPP_

#include "clODE_struct_defs.cl"
#include "OpenCLResource.hpp"

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include "OpenCL/cl2.hpp"

#include <map>
#include <string>
#include <vector>

struct ProblemInfo
{
    std::string clRHSfilename;
    cl_int nVar;
    cl_int nPar;
    cl_int nAux;
    cl_int nWiener;
    std::vector<std::string> varNames;
    std::vector<std::string> parNames;
    std::vector<std::string> auxNames;
};

class CLODE
{

protected:
    //Problem details
    ProblemInfo prob;
    std::string clRHSfilename;
    cl_int nVar, nPar, nAux, nWiener;
    cl_int nPts = 1;

    //Stepper specification
    std::string stepper;
    std::vector<std::string> availableSteppers;
    std::map<std::string, std::string> stepperDefineMap;

    bool clSinglePrecision;
    size_t realSize;

    //Compute device(s)
    OpenCLResource opencl;
    std::string clodeRoot;

    cl_int nRNGstate = 2; //TODO: different RNGs could be selected like steppers...?

    SolverParams<cl_double> sp;
    std::vector<cl_double> tspan, x0, pars, xf, dt;
    size_t x0elements, parselements, RNGelements;

    std::vector<cl_ulong> RNGstate;

    //Device variables
    cl::Buffer d_tspan, d_x0, d_pars, d_sp, d_xf, d_RNGstate, d_dt;

    //kernel object
    std::string clprogramstring, buildOptions, ODEsystemsource;
    cl::Kernel cl_transient;

    //flag to indicate whether kernel can be executed
    bool clInitialized = false;

    
    void setCLbuildOpts(std::string extraBuildOpts = "");
    std::string getStepperDefine();
    SolverParams<cl_float> solverParamsToFloat(SolverParams<cl_double> sp);

    //~private:
    //~ CLODE( const CLODE& other ); // non construction-copyable
    //~ CLODE& operator=( const CLODE& ); // non copyable


public:
    //for now, require all arguments. TODO: convenience constructors?
    CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot);
    CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot);
    virtual ~CLODE();

    //Set functions: trigger rebuild etc
    void setNewProblem(ProblemInfo prob);               //requires rebuild: pars/vars. Opencl context OK
    void setStepper(std::string newStepper);            //requires rebuild: Host + Device data OK
    void setPrecision(bool clSinglePrecision);          //requires rebuild: all device vars. Opencl context OK
    void setOpenCL(OpenCLResource opencl);              //requires rebuild: all device vars. Host problem data OK
    void setOpenCL(unsigned int platformID, unsigned int deviceID);
    void setClodeRoot(const std::string clodeRoot);

    void buildProgram(std::string extraBuildOpts = ""); //build the program object (inherited by subclasses)
    void buildCL(); // build program and create kernel objects - overloaded by subclasses to include any extra kernels

    // set all problem data needed to run
    virtual void initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp);

    void setNpts(cl_int newNpts); //resizes the nPts-dependent input variables
    void setProblemData(std::vector<cl_double> newX0, std::vector<cl_double> newPars); //set both pars and X0 to change nPts
    void setTspan(std::vector<cl_double> newTspan);
    void setX0(std::vector<cl_double> newX0);     //no change in nPts: newX0 must match nPts
    void setPars(std::vector<cl_double> newPars); //no change in nPts: newPars must match nPts
    void setSolverParams(SolverParams<cl_double> newSp);

    void seedRNG();
    void seedRNG(cl_int mySeed); //overload for setting reproducible seeds

    //simulation routine. TODO: overloads?
    void transient(); //integrate forward using stored tspan, x0, pars, and solver pars
    // void transient(std::vector<cl_double> newTspan); //integrate forward using stored x0, pars, and solver pars
    // void transient(std::vector<cl_double> newTspan, std::vector<cl_double> newX0); //integrate forward using stored pars, and solver pars
    // void transient(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars); //integrate forward using stored solver pars
    // void transient(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp);

    void shiftTspan(); //t0 <- tf, tf<-(tf + tf-t0)
    void shiftX0();    //d_x0 <- d_xf (device to device transfer)

    std::vector<cl_double> getTspan() { return tspan; };
    std::vector<cl_double> getX0();
    std::vector<cl_double> getXf();
    std::string getProgramString();
    std::vector<std::string> getAvailableSteppers() { return availableSteppers; };

    void printStatus();

    // Getters
    const ProblemInfo getProblemInfo() const { return prob; };
};

#endif //CLODE_HPP_
