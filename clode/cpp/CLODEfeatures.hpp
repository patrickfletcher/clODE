//
// Created by Patrick Fletcher 2017
//

#ifndef CLODE_FEATURES_HPP_
#define CLODE_FEATURES_HPP_

#include "CLODE.hpp"
#include "clODE_struct_defs.cl"
#include "observers.cl"
#include "OpenCLResource.hpp"

#define CL_HPP_ENABLE_EXCEPTIONS
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#define CL_HPP_ENABLE_PROGRAM_CONSTRUCTION_FROM_ARRAY_COMPATIBILITY
#include "OpenCL/cl2.hpp"

#include <map>
#include <string>
#include <vector>

class CLODEfeatures : public CLODE
{

protected:
    std::string observer;
    size_t ObserverParamsSize;

    std::map<std::string, struct ObserverInfo> observerDefineMap;
    std::vector<std::string> featureNames;
    std::vector<std::string> availableObserverNames;

    int nFeatures;
    size_t observerDataSize;
    std::vector<cl_double> F;
    ObserverParams<cl_double> op;
    size_t Felements;
    bool observerInitialized = false;

    cl::Buffer d_odata, d_op, d_F;
    cl::Kernel cl_initializeObserver;
    cl::Kernel cl_features;

    std::string observerBuildOpts;
    std::string observerName;

    ObserverParams<cl_float> observerParamsToFloat(ObserverParams<cl_double> op);
    std::string getObserverBuildOpts();
    void updateObserverDefineMap(); // update host variables representing feature detector: nFeatures, featureNames, observerDataSize
    void resizeFeaturesVariables(); // d_odata and d_F depend on nPts. nPts change invalidates d_odata

public:
    CLODEfeatures(ProblemInfo prob, std::string stepper, std::string observer, ObserverParams<cl_double> op, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot);
    CLODEfeatures(ProblemInfo prob, std::string stepper, std::string observer, ObserverParams<cl_double> op, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot);
    virtual ~CLODEfeatures();

    void buildCL() override; // build program and create kernel objects

    void setObserverParams(ObserverParams<cl_double> newOp); // requires rebuild: maxEventTimestamps in ObserverData... TODO: ideally decouple this!
    void setObserver(std::string newObserver);               // requires rebuild: program, kernel, kernel args. Host + Device data OK

    // simulation routine.
    void initializeObserver();                 // initialize Observer struct: possibly integrate forward an interval of duration (tf-t0), rewinds to t0
    void features();                           // integrate forward using stored tspan, x0, pars, and solver pars
    void features(bool reinitialize_observer); // allow manually forcing re-init of observer data
    bool isObserverInitialized() { return observerInitialized; };

    // Get functions
    const ObserverParams<cl_double> getObserverParams() const { return op; };
    const std::string getObserverName() const { return observerName; }
    const std::vector<cl_double> getF();
    const int getNFeatures() const { return nFeatures; };
    const std::vector<std::string> getFeatureNames() const { return featureNames; };
    const std::vector<std::string> getAvailableObservers() const { return availableObserverNames; };
};

#endif // CLODE_FEATURES_HPP_
