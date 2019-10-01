#include <stdexcept>
#include <stdio.h>
#include <algorithm>
#include <string>

#include "CLODEtrajectory.hpp"

CLODEtrajectory::CLODEtrajectory(ProblemInfo prob, StepperType stepper, bool clSinglePrecision, OpenCLResource opencl) 
    : CLODE(prob, stepper, clSinglePrecision, opencl), nStoreMax(0) {
		
    //printf("\nCLODE has been specialized: CLODEtrajectory\n");
	
	clprogramstring+=read_file(clodeRoot+"trajectory.cl");
}

CLODEtrajectory::~CLODEtrajectory(){}


//initialize everything
void CLODEtrajectory::initialize(std::vector<double> newTspan, std::vector<double> newX0, std::vector<double> newPars, SolverParams<double> newSp) {
	
	//(re)build the program
	buildProgram("");
	
	//Base class initialize kernel:
	initializeTransientKernel();
    
	//initialize the trajectory kernel
	initializeTrajectoryKernel();
    
	setSolverParams(newSp);
    setTspan(newTspan);
    
	setProblemData(newX0, newPars); //will set nPts
	
	resizeTrajectoryVariables(); //set up output variables too, which depend on nPts and nStoreMax=(tspan[1]-tspan[0])/(sp.dt*sp.nout)+1
	
	clInitialized=true;
}


void CLODEtrajectory::initializeTrajectoryKernel() {
	
	try { 
        cl_trajectory = cl::Kernel(opencl.getProgram(), "trajectory", &opencl.error);
    }
    catch (cl::Error &er) {
        printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
        throw er;
    }
}


void CLODEtrajectory::resizeTrajectoryVariables() {

	//catch any changes to tspan, dt, or nout:
	int currentStoreAlloc=(tspan[1]-tspan[0])/(sp.dt*sp.nout)+1;
	
	currentStoreAlloc=std::min(currentStoreAlloc, sp.max_store);
	//~ if (currentStoreAlloc > sp.max_store) {
		//~ currentStoreAlloc=sp.max_store;
		//~ printf("Warning: (tspan[1]-tspan[0])/(dt*nout)+1 > max_store. Reducing storage to max_store=%d timepoints. \n",sp.max_store);
	//~ }
	
	//check largest desired memory chunk against device's maximum allowable variable size
	size_t largestAlloc = std::max(1/*t*/, std::max(nVar/*x, dx*/, nAux/*aux*/))*nPts*currentStoreAlloc*realSize;
	
	if (largestAlloc > opencl.getMaxMemAllocSize()) {
		currentStoreAlloc=std::floor(opencl.getMaxMemAllocSize()/(std::max(1, std::max(nVar, nAux))*nPts*realSize));
		printf("ERROR: storage requested exceeds device maximum vairable size. Reason: %s. Try reducing storage to <=%d time points, or reducing nPts. \n",nAux>nVar?"aux vars":"state vars", currentStoreAlloc);
		throw std::invalid_argument("nPts*nStoreMax*nVar*realSize or nPts*nStoreMax*nAux*realSize is too big");
	}
	
	size_t currentTelements=currentStoreAlloc * nPts;
	
	//only resize device variables if size changed, or if not yet initialized
	if (!clInitialized || nStoreMax!=currentStoreAlloc || telements!=currentTelements) {
		
		nStoreMax=currentStoreAlloc;
		telements=currentTelements;
		xelements=nVar * currentTelements;
		auxelements=nAux * currentTelements;
		
		t.resize(telements);
		x.resize(xelements);
		dx.resize(xelements);
		aux.resize(auxelements);
		nStored.resize(nPts);
		
		//resize device variables
		try {
			//trajectory
			d_t  = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * telements, NULL, &opencl.error); 
			d_x  = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * xelements, NULL, &opencl.error); 
			d_dx  = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * xelements, NULL, &opencl.error); 
			d_aux  = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * auxelements, NULL, &opencl.error); 
			d_nStored = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, sizeof(int)*nPts, NULL, &opencl.error); 
			
		}
		catch (cl::Error &er) {
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
			throw er;
		}
	}
}



//Simulation routine
void CLODEtrajectory::trajectory() {
    if (clInitialized) {
		//resize output variables - will only occur if nPts or nSteps has changed [~4ms overhead on Tornado]
		resizeTrajectoryVariables();
			
		try {
			//kernel arguments
			cl_trajectory.setArg(0, d_tspan);
			cl_trajectory.setArg(1, d_x0); 
			cl_trajectory.setArg(2, d_pars); 
			cl_trajectory.setArg(3, d_sp);
			cl_trajectory.setArg(4, d_xf); 
			cl_trajectory.setArg(5, d_RNGstate); 
			cl_trajectory.setArg(6, d_t); 
			cl_trajectory.setArg(7, d_x); 
			cl_trajectory.setArg(8, d_dx); 
			cl_trajectory.setArg(9, d_aux); 
			cl_trajectory.setArg(10, nStoreMax);
			cl_trajectory.setArg(11, d_nStored); 
			
			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_trajectory, cl::NullRange, cl::NDRange(nPts), cl::NullRange);
			opencl.getQueue().finish();
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



std::vector<double> CLODEtrajectory::getT() {

		if (clSinglePrecision) { //cast back to double
			std::vector<float> tF(telements);
			opencl.error = copy(opencl.getQueue(), d_t, tF.begin(), tF.end());
			t.assign(tF.begin(),tF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_t, t.begin(), t.end());
		}
		
		return t;
}


std::vector<double> CLODEtrajectory::getX() {

		if (clSinglePrecision) { //cast back to double
			std::vector<float> xF(xelements);
			opencl.error = copy(opencl.getQueue(), d_x, xF.begin(), xF.end());
			x.assign(xF.begin(),xF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_x, x.begin(), x.end());
		}
		
		return x;
}


std::vector<double> CLODEtrajectory::getDx() {

		if (clSinglePrecision) { //cast back to double
			std::vector<float> dxF(xelements);
			opencl.error = copy(opencl.getQueue(), d_dx, dxF.begin(), dxF.end());
			dx.assign(dxF.begin(),dxF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_dx, dx.begin(), dx.end());
		}
		
		return dx;
}


std::vector<double> CLODEtrajectory::getAux() { //cast back to double

		if (clSinglePrecision) { 
			std::vector<float> auxF(auxelements);
			opencl.error = copy(opencl.getQueue(), d_aux, auxF.begin(), auxF.end());
			aux.assign(auxF.begin(),auxF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_aux, aux.begin(), aux.end());
		}
		
		return aux;
}

std::vector<int> CLODEtrajectory::getNstored() { 

		opencl.error = copy(opencl.getQueue(), d_nStored, nStored.begin(), nStored.end());
		return nStored;
}
