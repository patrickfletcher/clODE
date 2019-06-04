/* clODE: a simulator class to run parallel ODE simulations on OpenCL capable hardware.
 * A clODE simulator solves on initial value problem over a grid of parameters and/or initial conditions. At each timestep,
 * "observer" rountine may be called to record/store/compute features of the solutions. Examples include storing the full
 * trajectory, recording the times and values of local extrema in a variable of the system, or directly computing other 
 * features of the trajectory.  
 */

//TODO: namespaces?

//TODO: break up computation (timespan) into chunks that don't crash the system...

//TODO: separate build/initialize kernel from initialize?

//TODO: choosing specific RNG - get nRNGstate using a switch

//TODO: device-to-device transfers instead of overwriting x0?

//when compiling, be sure to provide the clODE root directory as a define:
// -DCLODE_ROOT="path/to/my/clODE/" 

#ifndef CLODE_HPP_
#define CLODE_HPP_

#define dbg_printf printf
//~ #define dbg_printf 

//if we are compiling from matlab MEX, redefine printf to mexPrintf so it prints to matlab command window.
#ifdef MATLAB_MEX_FILE
    #include "mex.h"
    #define printf mexPrintf
#endif

#include "OpenCLResource.hpp"
#include "clODE_struct_defs.h"


// some utility structs to pass to the OpenCL kernels.  TEMPLATE for float/double after debugging of pure double case.
struct ProblemInfo{
    std::string clRHSfilename;
    int nVar;
    int nPar;
    int nAux;
    int nWiener;
};


//enum for stepper selection
typedef enum StepperType{
	euler=0, 
	heun, 
	rungeKutta4,
	heunEuler,
	bogackiShampine23,
	dorpri5,
} StepperType;



class CLODE {
	
	protected:
	
		//Problem details (from ProblemInfo struct)
        std::string clRHSfilename; 
		int nVar, nPar, nAux, nWiener; 
		cl_int nPts;
		
		//Stepper specification
        StepperType stepper;
        std::string stepperName;
        
        bool clSinglePrecision;
        size_t realSize;
        
        //Compute device(s)
        OpenCLResource opencl;
        const std::string clodeRoot;
		
        int nRNGstate; //TODO: different RNGs could be selected like steppers...?
        
        SolverParams<double> sp;
        std::vector<double> tspan, x0, pars, auxf;
        size_t x0elements, parselements, auxfelements, RNGelements;
        
        std::vector<cl_ulong> RNGstate;
        
        //Device variables
        cl::Buffer d_tspan, d_x0, d_pars, d_sp, d_auxf, d_RNGstate; 

		//kernel object
        std::string clprogramstring, buildOptions, ODEsystemsource;
		cl::Kernel cl_transient;
        
        //flag system to ensure kernel can be executed
        bool clInitialized;


		std::string getStepperDefine();
		SolverParams<float> solverParamsToFloat(SolverParams<double> sp);
		
	//~private:
		//~ CLODE( const CLODE& other ); // non construction-copyable
		//~ CLODE& operator=( const CLODE& ); // non copyable
		
		void initializeTransientKernel();
        void setNpts(int newNpts);  //resizes the nPts-dependent input variables 
		
    public:
        
        //for now, require all arguments. TODO: convenience constructors?
        CLODE(ProblemInfo prob, StepperType stepper, bool clSinglePrecision, OpenCLResource opencl);
        //~ CLODE(); //must follow with all set functions to use
        //~ CLODE(ProblemInfo prob); //set stepper, precision, and opencl
        //~ CLODE(ProblemInfo prob, StepperType stepper=rungeKutta4, bool clSinglePrecision=true, OpenCLResource opencl=OpenCLResource()); //alt: use defaults?
        ~CLODE();
        
        //Set functions: trigger rebuild etc
		void setNewProblem(ProblemInfo prob);		//rebuild: program, kernel, kernel args, device vars. Opencl context OK
		void setStepper(StepperType newStepper);	//rebuild: program, kernel, kernel args. Host + Device data OK
		void setPrecision(bool clSinglePrecision);	//rebuild: program, kernel, kernel args, device vars. Opencl context OK
		void setOpenCL(OpenCLResource opencl);		//rebuild: context, program, kernel, kernel args, device vars. Host problem data OK
		void buildProgram(std::string extraBuildOpts="");      //called by initialize, to all change in prob/stepper/precision/openclResource

		//build program, set all problem data needed to run
        virtual void initialize(std::vector<double> newTspan, std::vector<double> newX0, std::vector<double> newPars, SolverParams<double> newSp);
        
        void setProblemData(std::vector<double> newX0, std::vector<double> newPars); //this routine may change nPts
        void setTspan(std::vector<double>  newTspan);
        void setX0(std::vector<double>  newX0); //no change in nPts
        void setPars(std::vector<double> newPars); //no change in nPts
        void setSolverParams(SolverParams<double> newSp);
        
		void seedRNG();
		void seedRNG(int mySeed); //TODO: overload for setting reproducible seeds
        
        //simulation routine
        void transient();  //integrate forward an interval of duration (tf-t0). automatically sets device tspan and x0 to be ready to continue forward from (tf,xf)
        
        std::vector<double> getTspan(){return tspan;};  //returns "current" tspan
        std::vector<double> getX0(); //returns the vector of x values at current time (x0 on device is overwritten with xf after integration)
        std::vector<double> getAuxf();
        std::string programString() {return clprogramstring;};
        
        void printStatus();
       
};


#endif //CLODE_HPP_
