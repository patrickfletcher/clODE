//
// Created by Patrick Fletcher 2017
//

#ifndef CLODE_HPP_
#define CLODE_HPP_

#include "clODE_struct_defs.cl"
#include "OpenCLResource.hpp"

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include "OpenCL/cl2.hpp"
// #include "OpenCL/opencl.hpp"

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

    ProblemInfo(std::string clRHSfilename, cl_int nVar, cl_int nPar, cl_int nAux, cl_int nWiener, std::vector<std::string> varNames, std::vector<std::string> parNames, std::vector<std::string> auxNames)
    {
        this->clRHSfilename = clRHSfilename;
        this->nVar = nVar;
        this->nPar = nPar;
        this->nAux = nAux;
        this->nWiener = nWiener;
        this->varNames = varNames;
        this->parNames = parNames;
        this->auxNames = auxNames;
    }

    ProblemInfo(std::string clRHSfilename, std::vector<std::string> varNames, std::vector<std::string> parNames, std::vector<std::string> auxNames = std::vector<std::string>(), cl_int nWiener = 0)
    {
        this->clRHSfilename = clRHSfilename;
        this->nVar = varNames.size();
        this->nPar = parNames.size();
        this->nAux = auxNames.size();
        this->nWiener = nWiener;
        this->varNames = varNames;
        this->parNames = parNames;
        this->auxNames = auxNames;
    }

    ProblemInfo()
    {
        this->clRHSfilename = "";
        this->nVar = 0;
        this->nPar = 0;
        this->nAux = 0;
        this->nWiener = 0;
        this->varNames = std::vector<std::string>();
        this->parNames = std::vector<std::string>();
        this->auxNames = std::vector<std::string>();
    }

    // These setters and getters exist for compatibility with the Python interface
    void setVarNames(std::vector<std::string> varNamesIn)
    {
        this->varNames = varNamesIn;
        this->nVar = varNamesIn.size();
    }

    void setParNames(std::vector<std::string> parNamesIn)
    {
        this->parNames = parNamesIn;
        this->nPar = parNamesIn.size();
    }

    void setAuxNames(std::vector<std::string> auxNamesIn)
    {
        this->auxNames = auxNamesIn;
        this->nAux = auxNamesIn.size();
    }

    std::vector<std::string> getVarNames() { return this->varNames; }
    std::vector<std::string> getParNames() { return this->parNames; }
    std::vector<std::string> getAuxNames() { return this->auxNames; }
};

class CLODE
{

protected:
    // Problem details
    ProblemInfo prob;
    std::string clRHSfilename;
    cl_int nVar, nPar, nAux, nWiener;
    cl_int nPts = 0;

    // Stepper specification
    std::string stepper;
    std::vector<std::string> availableSteppers;
    std::map<std::string, std::string> stepperDefineMap;

    bool clSinglePrecision;
    size_t realSize;

    // Compute device(s)
    OpenCLResource opencl;
    std::string clodeRoot;

    cl_int nRNGstate = 2; // TODO: different RNGs could be selected like steppers...?

    SolverParams<cl_double> sp;
    std::vector<cl_double> tspan, x0, pars, xf, dt;
    size_t x0elements, parselements, RNGelements;

    std::vector<cl_ulong> RNGstate;

    // Device variables
    cl::Buffer d_tspan, d_sp, d_pars, d_x0, d_xf, d_RNGstate, d_dt;

    // kernel object
    std::string clprogramstring, buildOptions, ODEsystemsource;
    cl::Kernel cl_transient;

    void setClodeRoot(const std::string clodeRoot);
    void setCLbuildOpts(std::string extraBuildOpts = "");
    void buildProgram(std::string extraBuildOpts = ""); // build the program object (inherited by subclasses)
    void setNpts(cl_int newNpts);                       // resizes the nPts-dependent input variables - can't be called alone as it doesn't populate the arrays with data
    std::string getStepperDefine();
    SolverParams<cl_float> solverParamsToFloat(SolverParams<cl_double> sp);

    //~private:
    //~ CLODE( const CLODE& other ); // non construction-copyable
    //~ CLODE& operator=( const CLODE& ); // non copyable

public:
    // for now, require all arguments. TODO: convenience constructors?
    CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot);
    CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot);
    virtual ~CLODE();

    // Setters that trigger need to rebuild CL program
    void setProblemInfo(ProblemInfo prob);     // Opencl context OK
    void setStepper(std::string newStepper);   // Host + Device data OK
    void setPrecision(bool clSinglePrecision); // reset all device vars. Opencl context OK
    void setOpenCL(OpenCLResource opencl);     // reset all device vars. Host problem data OK
    void setOpenCL(unsigned int platformID, unsigned int deviceID);

    virtual void buildCL(); // build program and create kernel objects - overloaded by subclasses to include any extra kernels

    // no need to rebuild CL program
    void setProblemData(std::vector<cl_double> newX0, std::vector<cl_double> newPars); // set both pars and X0 to change nPts
    void setTspan(std::vector<cl_double> newTspan);
    void setX0(std::vector<cl_double> newX0);     // no change in nPts: newX0 must match nPts
    void setPars(std::vector<cl_double> newPars); // no change in nPts: newPars must match nPts
    void setSolverParams(SolverParams<cl_double> newSp);

    void seedRNG();
    void seedRNG(cl_int mySeed); // overload for setting reproducible seeds

    // simulation routine. TODO: overloads?
    void transient(); // integrate forward using stored tspan, x0, pars, and solver pars
    // void transient(std::vector<cl_double> newTspan); //integrate forward using stored x0, pars, and solver pars
    // void transient(std::vector<cl_double> newTspan, std::vector<cl_double> newX0); //integrate forward using stored pars, and solver pars
    // void transient(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars); //integrate forward using stored solver pars
    // void transient(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp);

    void shiftTspan(); // t0 <- tf, tf<-(tf + tf-t0)
    void shiftX0();    // d_x0 <- d_xf (device to device transfer)

    // Getters
    const ProblemInfo getProblemInfo() const { return prob; };
    const std::vector<cl_double> getTspan() const { return tspan; };
    const SolverParams<cl_double> getSolverParams() const { return sp; };
    const std::vector<cl_double> getPars() const { return pars; };
    const std::vector<cl_double> getX0();
    const std::vector<cl_double> getXf();
    const std::vector<std::string> getAvailableSteppers() const { return availableSteppers; };

    const std::string getProgramString() const { return buildOptions + clprogramstring + ODEsystemsource; };
    void printStatus();
};

#endif // CLODE_HPP_
