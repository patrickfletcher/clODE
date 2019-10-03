#define dbg_printf printf
//~ #define dbg_printf 

//if we are compiling from matlab MEX, redefine printf to mexPrintf so it prints to matlab command window.
#ifdef MATLAB_MEX_FILE
    #include "mex.h"
    #define printf mexPrintf
#endif

#include <vector>
#include <stdio.h>
#include <string>
#include <algorithm> //std::max

#include "CLODEfeatures.hpp"

CLODEfeatures::CLODEfeatures(ProblemInfo prob, StepperType stepper, ObserverType observer, bool clSinglePrecision, OpenCLResource opencl) 
    : CLODE(prob, stepper, clSinglePrecision, opencl), observer(observer)
{
		
    //printf("\nCLODE has been specialized: CLODEfeatures\n");
		
	clprogramstring+=read_file(clodeRoot+"features.cl");
}

CLODEfeatures::~CLODEfeatures(){}


//initialize everything
void CLODEfeatures::initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp, ObserverParams<cl_double> newOp) {

	clInitialized = false;
	observerBuildOpts=getObserverBuildOpts();
	
	//(re)build the program
	buildProgram(observerBuildOpts);
    // printf("Using observer: %s\n",observerName.c_str());
	
	//Base class initialize kernel:
	initializeTransientKernel();
	// printf("Transient kernel initialized.\n");
    
	//initialize the trajectory kernel
	initializeFeaturesKernel();
	// printf("Features kernel initialized.\n");
    
	setSolverParams(newSp);
	// printf("Solver parameters set.\n");
	setObserverParams(newOp);
	// printf("Observer parameters set.\n");
    setTspan(newTspan);
	// printf("Integration time span set.\n");
    
	setProblemData(newX0, newPars); //will set nPts
	// printf("Problem data set.\n");
	
	resizeFeaturesVariables(); //set up d_F and d_odata too, which depend on nPts
	
	clInitialized=true;
	// printf("Initialization complete.\n");
	
}


std::string CLODEfeatures::getObserverBuildOpts() {
	std::string observerdefine="";
	switch (observer) {
        case basic:
            observerdefine=" -DOBSERVER_BASIC ";
            observerName="OBSERVER_BASIC";
            nFeatures=6;
            observerDataSize=5*realSize+2*sizeof(cl_int);
            observerDataSize=observerDataSize+observerDataSize%realSize; //need to align the struct to a multiple of the largest type inside (here, realtype)...
            break;
            
        case basicAllVar:
            observerdefine=" -DOBSERVER_BASIC_ALLVAR ";
            observerName="OBSERVER_BASIC_ALLVAR";
            nFeatures=5*nVar+3*nAux+1;
            observerDataSize=5*realSize*nVar+3*realSize*nAux+2*sizeof(cl_int);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
            break;
            
        case localmax:
            observerdefine=" -DOBSERVER_LOCAL_MAX ";
            observerName="OBSERVER_LOCAL_MAX";
            nFeatures=18;
            observerDataSize=29*realSize+3*sizeof(cl_int);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
            break;
            
        case section1:
            observerdefine=" -DOBSERVER_SECTION_1 ";
            observerName="OBSERVER_SECTION_1";
            nFeatures=30;
            observerDataSize=2*sizeof(cl_int);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
            break;
            
        case neighborhood1:
            observerdefine=" -DOBSERVER_NEIGHBORHOOD_1 ";
            observerName="OBSERVER_NEIGHBORHOOD_1";
            nFeatures=30;
            observerDataSize=2*sizeof(cl_int);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
            break;
            
        case section2:
            observerdefine=" -DOBSERVER_SECTION_2 ";
            observerName="OBSERVER_SECTION_2";
            nFeatures=11;
            observerDataSize=(20+6*nVar)*realSize+4*sizeof(cl_int)+sizeof(cl_bool);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
            break;
            
        case neighborhood2:
            observerdefine=" -DOBSERVER_NEIGHBORHOOD_2 ";
            observerName="OBSERVER_NEIGHBORHOOD_2";
            nFeatures=11;
            observerDataSize=(21+7*nVar)*realSize+3*sizeof(cl_int)+3*sizeof(cl_bool);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
			printf("observerDataSize = %d", observerDataSize);
            break; // ObserverData size: 3*int + (4*3 + 7*nVar + 9)*realtype + 3*bool
            
        default:
            observerdefine=" -DOBSERVER_BASIC ";
            observerName="OBSERVER_BASIC";
            nFeatures=6;
            observerDataSize=5*realSize+2*sizeof(cl_int);
            observerDataSize=observerDataSize+observerDataSize%realSize; 
            break;
    }
    return observerdefine;
}


void CLODEfeatures::setObserver(ObserverType newObserver) {
	if (newObserver!=observer) {
		observer=newObserver;
		clInitialized=false;
	}
}


void CLODEfeatures::initializeFeaturesKernel() {
	
	try { 
		if (clSinglePrecision) { //downcast to float if desired
			d_op = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(ObserverParams<cl_float>), NULL, &opencl.error); 
		}
		else {
			d_op = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(ObserverParams<cl_double>), NULL, &opencl.error); 
		}
        cl_features = cl::Kernel(opencl.getProgram(), "features", &opencl.error);
        
		size_t preferred_multiple;
		cl::Device dev;
		opencl.getProgram().getInfo(CL_PROGRAM_DEVICES,&dev);
		cl_features.getWorkGroupInfo(dev,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,&preferred_multiple);
		printf("Preferred work size multiple (features): %d\n",preferred_multiple);
    }
    catch (cl::Error &er) {
        printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
        throw er;
    }
}


void CLODEfeatures::setObserverParams(ObserverParams<cl_double> newOp){
    
    op=newOp;
    try {
		if (clSinglePrecision) { //downcast to float if desired
			ObserverParams<cl_float> opF=observerParamsToFloat(op);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_op, CL_TRUE, 0, sizeof(opF), &opF);
		}
		else {
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_op, CL_TRUE, 0, sizeof(op), &op);
		}
    }
    catch (cl::Error &er) {
        printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
        throw er;
    }
}


//TODO: define an assignment/type cast operator in the struct?
ObserverParams<cl_float> CLODEfeatures::observerParamsToFloat(ObserverParams<cl_double> op) {
	ObserverParams<cl_float> opF;
	
	opF.eVarIx=op.eVarIx;
	opF.fVarIx=op.fVarIx;
	opF.maxEventCount=op.maxEventCount;
	opF.minXamp=op.minXamp;
	opF.nHoodRadius=op.nHoodRadius;
	opF.xUpThresh=op.xUpThresh;
	opF.xDownThresh=op.xDownThresh;
	opF.dxUpThresh=op.dxUpThresh;
	opF.dxDownThresh=op.dxDownThresh;
	opF.eps_dx=op.eps_dx;
	
	return opF;
}




void CLODEfeatures::resizeFeaturesVariables() {

	size_t currentFelements=nFeatures * nPts;
	
	size_t largestAlloc = std::max(nFeatures*realSize, observerDataSize) * nPts;
	
	if (largestAlloc > opencl.getMaxMemAllocSize()) {
		int maxNpts=floor(opencl.getMaxMemAllocSize()/(cl_ulong)std::max(nFeatures*realSize, observerDataSize));
		printf("nPts is too large, requested memory size exceeds selected device's limit. Maximum nPts appears to be %d \n", maxNpts);
		throw std::invalid_argument("nPts is too large");
	}
	
	//resize device variables if nPts changed, or if not yet initialized
	if (!clInitialized || Felements!=currentFelements) {
		
		Felements=currentFelements;
		F.resize(currentFelements);
		
		//resize device variables
		try {
			d_odata  = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, observerDataSize * nPts, NULL, &opencl.error);
			d_F  = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * currentFelements, NULL, &opencl.error); 
		}
		catch (cl::Error &er) {
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
			throw er;
		}
	}
}


//Simulation routines

//overload to allow manual re-initialization of observer data at any time.
void CLODEfeatures::features(bool newDoObserverInitFlag) {
	doObserverInitialization=newDoObserverInitFlag;
	features();
}


//
void CLODEfeatures::features() {
	
    if (clInitialized) {
		
		//resize output variables - will only occur if nPts has changed 
		resizeFeaturesVariables();
			
		try {
			//kernel arguments
			cl_features.setArg(0, d_tspan);
			cl_features.setArg(1, d_x0); 
			cl_features.setArg(2, d_pars); 
			cl_features.setArg(3, d_sp);
			cl_features.setArg(4, d_xf); 
			cl_features.setArg(5, d_RNGstate); 
			cl_features.setArg(6, d_odata); 
			cl_features.setArg(7, d_op); 
			cl_features.setArg(8, d_F); 
			cl_features.setArg(9, doObserverInitialization);
			
			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_features, cl::NullRange, cl::NDRange(nPts), cl::NullRange);
			opencl.getQueue().finish();
			
			doObserverInitialization=0;
		}
		catch (cl::Error &er) {
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
			throw er;
		}
	}
	else {
		printf("CLODE has not been initialized\n");
	}
}



std::vector<double> CLODEfeatures::getF() {

		if (clSinglePrecision) { //cast back to double
			std::vector<cl_float> FF(Felements);
			opencl.error = copy(opencl.getQueue(), d_F, FF.begin(), FF.end());
			F.assign(FF.begin(),FF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_F, F.begin(), F.end());
		}
		
		return F;
}

