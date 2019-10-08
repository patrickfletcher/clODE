/* clODE: a simulator class to run parallel ODE simulations on OpenCL capable hardware.
 * A clODE simulator solves on initial value problem over a grid of parameters and/or initial conditions. At each timestep,
 * "observer" rountine may be called to record/store/compute features of the solutions. Examples include storing the full
 * trajectory, recording the times and values of local extrema in a variable of the system, or directly computing other 
 * features of the trajectory.  
 */


//when compiling, be sure to provide the clODE root directory as a define:
// -DCLODE_ROOT="path/to/my/clODE/" 

#ifndef CLODE_FEATURES_HPP_
#define CLODE_FEATURES_HPP_

#include "CLODE.hpp"
#include "clODE_struct_defs.h"
#include "observers.h"
#include "OpenCLResource.hpp"

#define __CL_ENABLE_EXCEPTIONS
#if defined(__APPLE__) || defined(__MACOSX)
    #include "OpenCL/cl.hpp"
#else
    #include <CL/cl.hpp>
#endif

#include <string>
#include <vector>


//enum for observer type
typedef enum ObserverType{
	basic=0,
	basicAllVar,
	localmax,
	section1,
	neighborhood1, 
	section2,
	neighborhood2
} ObserverType;


class CLODEfeatures: public CLODE {
		
	protected:
		
		ObserverType observer;
		size_t ObserverParamsSize;
		
        int nFeatures;
        size_t observerDataSize;
        std::vector<cl_double> F;
        ObserverParams<cl_double> op;
        size_t Felements;
		int doObserverInitialization=1;
		
        cl::Buffer d_odata, d_op, d_F; 
		cl::Kernel cl_features;
		
		std::string observerBuildOpts;
		std::string observerName;
		
		ObserverParams<float> observerParamsToFloat(ObserverParams<double> sp);
		
		std::string getObserverBuildOpts();
		void initializeFeaturesKernel();
		void resizeFeaturesVariables(); //d_odata and d_F depend on nPts. nPts change invalidates d_odata
		
    public:
        
        CLODEfeatures(ProblemInfo prob, StepperType stepper, ObserverType observer, bool clSinglePrecision, OpenCLResource opencl);
        ~CLODEfeatures();

		//build program, set all problem data needed to run
        virtual void initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp, ObserverParams<cl_double> newOp);
        
        void setObserverParams(ObserverParams<cl_double> newOp);
		void setObserver(ObserverType newObserver);	//rebuild: program, kernel, kernel args. Host + Device data OK
        
        //simulation routines.
        //TODO: overload transient so that it sets doObserverInitialization=1
        //TODO: overload with newX0, newPars; all four? 
        //TODO: pre-features? update edata basics (min/max x dx) but no event function
        void features(bool newDoObserverInitFlag); //allow manually forcing re-init of observer data
        void features();  //integrate forward an interval of duration (tf-t0)
        
        //Get functions
        int getNFeatures() {return nFeatures;};
        std::vector<double> getF();   
};


#endif //CLODE_FEATURES_HPP_
