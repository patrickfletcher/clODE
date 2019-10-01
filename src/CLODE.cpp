#include <stdexcept>
#include <stdio.h>
#include <string>
#include <algorithm>
#include <random>
#include <stdint.h>

#include "OpenCLResource.hpp"  
#include "CLODE.hpp"
#include "clODE_struct_defs.h"

//constructor sets problem info and builds the base clprogramstring
CLODE::CLODE(ProblemInfo prob, StepperType stepper, bool clSinglePrecision, OpenCLResource opencl) 
    : clRHSfilename(prob.clRHSfilename),nVar(prob.nVar),nPar(prob.nPar),nAux(prob.nAux),nWiener(prob.nWiener), nPts(0),
    stepper(stepper), clSinglePrecision(clSinglePrecision), realSize(clSinglePrecision?sizeof(cl_float):sizeof(cl_double)),
    opencl(opencl), clodeRoot(CLODE_ROOT), nRNGstate(2)
{
    //printf("\nCLODE base class constructor\n");
	
    clprogramstring=read_file(clodeRoot+"transient.cl");
}


CLODE::~CLODE()
{
}


void CLODE::setNewProblem(ProblemInfo newProb) {
	
	clRHSfilename=newProb.clRHSfilename;
	nVar=newProb.nVar;
	nPar=newProb.nPar;
	nAux=newProb.nAux;
	nWiener=newProb.nWiener;
	
	clInitialized=false;
}


void CLODE::setStepper(StepperType newStepper) {
	if (newStepper!=stepper) {
		stepper=newStepper;
		clInitialized=false;
	}
}


void CLODE::setPrecision(bool newPrecision) {
	if (newPrecision!=clSinglePrecision) {
		clSinglePrecision=newPrecision;
		realSize=newPrecision?sizeof(cl_float):sizeof(cl_double);
		clInitialized=false;
	}
}


void CLODE::setOpenCL(OpenCLResource newOpencl) {
	//~ if (newOpencl!=opencl) {
		opencl=newOpencl;
		clInitialized=false;
	//~ }
}

//build creates build option defined constants based on selected options, adds the ODEsystem source to clprogramstring then builds for selected OpenCL resource
void CLODE::buildProgram(std::string extraBuildOpts) {

    
	if (!clSinglePrecision && !opencl.getDoubleSupport()) {
		clSinglePrecision=true;
		printf("Warning: device selected does not support double precision. Using single precision\n");
	}
   
    buildOptions="";
    
   //specify precision
    if (clSinglePrecision)
        buildOptions+=" -DCLODE_SINGLE_PRECISION";
    else
        buildOptions+=" -DCLODE_DOUBLE_PRECISION";
    
    //specify stepper
    buildOptions+=getStepperDefine();
    
    //specify problem dimensions
    buildOptions+=" -DN_PAR=" + std::to_string((long long)nPar); //for older c++ compilers the to_string(int) overload of the STL isn't present
	buildOptions+=" -DN_VAR=" + std::to_string((long long)nVar);
	buildOptions+=" -DN_AUX=" + std::to_string((long long)nAux);
	buildOptions+=" -DN_WIENER=" + std::to_string((long long)nWiener);
    
    //include folder for CLODE
    buildOptions+=" -I" + clodeRoot;
    
    buildOptions+=extraBuildOpts;
    
    //ODEsystem source is delayed to here, in case we change it
    ODEsystemsource=read_file(clRHSfilename);
    
    // printf("%s", clprogramstring.c_str());
    // printf("%s", ODEsystemsource.c_str());
    // printf("%s", buildOptions.c_str());
	
    //now build
    opencl.buildProgramFromString(clprogramstring+ODEsystemsource,buildOptions);
    
	
    printStatus();
}

    
//TODO make stepperDefine a map? 
std::string CLODE::getStepperDefine() {
	
	std::string stepperDefine;
	
    switch (stepper) {
        case euler:
            stepperDefine=" -DFORWARD_EULER";
            stepperName="Forward Euler";
            break;
        case heun:
            stepperDefine=" -DHEUN";
            stepperName="Heun's method";
            break;
        case rungeKutta4:
            stepperDefine=" -DRUNGE_KUTTA4";
            stepperName="Runge-Kutta 4";
            break;
        case heunEuler:
            stepperDefine=" -DHEUN_EULER";
            stepperName="Adaptive Heun-Euler";
            break;
        case bogackiShampine23:
            stepperDefine=" -DBOGACKI_SHAMPINE23";
            stepperName="Bogacki-Shampine 23";
            break;
        case dorpri5:
            stepperDefine=" -DDORPRI5";
            stepperName="Dormand-Prince 45";
            break;
        default:
            stepperDefine=" -DRUNGE_KUTTA4";
            stepperName="Runge-Kutta 4";
    }
    
    return stepperDefine;
}


//initialize everything: build the program, create the kernels, and set all needed problem data.
void CLODE::initialize(std::vector<double> newTspan, std::vector<double> newX0, std::vector<double> newPars, SolverParams<double> newSp) {
	
	//(re)build the program
	buildProgram();
	
	//set up the kernel
    initializeTransientKernel();
    
	setSolverParams(newSp);
    setTspan(newTspan);
	setProblemData(newX0, newPars); //will call setNpts
	
	clInitialized=true;
}


void CLODE::initializeTransientKernel() {
	
	try { //declare device arrays that won't change size: tspan, SolverParams
		d_tspan = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, realSize*2, NULL, &opencl.error); 
		if (clSinglePrecision) {
			d_sp = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(SolverParams<cl_float>), NULL, &opencl.error); 
		}
		else
			d_sp = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(SolverParams<cl_double>), NULL, &opencl.error); 
			
		//initialize kernel and assign kernel arguments
        cl_transient = cl::Kernel(opencl.getProgram(), "transient", &opencl.error);
    }
    catch (cl::Error &er) {
        printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
        throw er;
    }
}


//initialize new set of trajectories (nPts may change)
void CLODE::setProblemData(std::vector<double> newX0, std::vector<double> newPars) {

    
	//check if newX0 and newPars are valid, and update nPts if needed:
	if (newX0.size()%nVar!=0) {
        printf("Invalid initial condition vector: not a multiple of nVar=%d\n", nVar);
        printf("...Initial conditions were not updated!\n");
        return;
	}
	
	if (newPars.size()%nPar!=0) {
        printf("Invalid parameter vector: not a multiple of nPar=%d\n", nPar);
        printf("...Parameters were not updated!\n");
        return;
	}

	// now check if newX0 and newPars represent same number of sets
	int nPtsX0=newX0.size()/nVar;
	int nPtsPars=newPars.size()/nPar;
	printf("Computed nPts: %d %d\n", nPtsX0, nPtsPars);
	if (nPtsX0!=nPtsPars) {
        printf("Initial contition and parameter vector dimensions don't match");
        printf("...Expected %d sets of each, recieved %d for x0 and %d for pars\n", nPts, nPtsX0, nPtsPars);
        printf("...Problem data was not updated!\n");
        return;
	}

	//set nPts
	setNpts(nPtsX0);
	printf("set nPts\n");

	//set things that depend on nPts
	setX0(newX0);
	printf("set X0\n");
	setPars(newPars);
	printf("set P\n");
}


//resize all the nPts dependent variables, only if nPts changed
void CLODE::setNpts(int newNpts){
	
	//unlikely that any of these should ever exceed memory limits...
	size_t largestAlloc = std::max(nVar, std::max(nPar, nAux)) * nPts * realSize;
	printf("Computed largestAlloc: %d\n", largestAlloc);
	
	if (largestAlloc > opencl.getMaxMemAllocSize()) {
		throw std::invalid_argument("nPts*nVar, nPts*nPar, or nPts*nAux is too large");
	}
	
	if (!clInitialized || newNpts!=nPts) {
		nPts=newNpts;
		
		x0elements=nVar * nPts;
		parselements=nPar * nPts;
		auxfelements=nAux * nPts; //TODO: handle nAux=0 case
		RNGelements=nRNGstate * nPts;

		//resize host variables
		x0.resize(x0elements);
		pars.resize(parselements);
		xf.resize(x0elements);
		auxf.resize(auxfelements);
		RNGstate.resize(RNGelements); 
		

		//resize device variables
		try {
			d_x0   = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * x0elements, NULL, &opencl.error);
			// printf("1\n");
			d_pars = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, realSize * parselements, NULL, &opencl.error);
			// printf("2\n");
			d_xf   = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * x0elements, NULL, &opencl.error);
			// printf("1\n");
			d_auxf  = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * auxfelements, NULL, &opencl.error);
			// printf("3\n");
			d_RNGstate  = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, sizeof(cl_ulong) * RNGelements, NULL, &opencl.error);
			// printf("4\n");
		}
		catch (cl::Error &er) {
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
			throw er;
		}
		
		//seed RNG must occur after device variable d_RNGstate is resized
		seedRNG();  
	}
}


void CLODE::setTspan( std::vector<double>  newTspan) {
	
	tspan=newTspan;
	
	//sync to device
	try {
		if (clSinglePrecision) { //downcast to float if desired
			std::vector<float> tspanF(tspan.begin(),tspan.end());
			opencl.error = copy(opencl.getQueue(), tspanF.begin(), tspanF.end(), d_tspan);
		}
		else {
			opencl.error = copy(opencl.getQueue(), tspan.begin(), tspan.end(), d_tspan);
		}
	}
	catch (cl::Error &er) {
		printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
		throw er;
	}
}

void CLODE::updateTspan(){
	std::vector<double> newTspan({tspan[1], tspan[1]+(tspan[1]-tspan[0])});
	setTspan(newTspan);
}


//set new x0. Cannot update nPts
void CLODE::setX0( std::vector<double>  newX0){
    
    if (newX0.size()==(size_t)nPts*nVar) {
        x0=newX0;
        
        //sync to device
        try {
			if (clSinglePrecision) { //downcast to float if desired
				std::vector<float> x0F(x0.begin(),x0.end());
				opencl.error = copy(opencl.getQueue(), x0F.begin(), x0F.end(), d_x0);
			}
			else {
				opencl.error = copy(opencl.getQueue(), x0.begin(), x0.end(), d_x0);
			}
        	
        }
        catch (cl::Error &er) {
            printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
            throw er;
        }        
    }
    else {
        printf("Invalid initial condition vector: Expected %d*%d elements, recieved %lu\n", nPts, nVar, newX0.size());
        printf("...Initial conditions were not updated!\n");
        //~ throw std::invalid_argument("Initial Condition vector has incorrect size.");
    }
}

void CLODE::updateX0(){
	//device to device transfer of Xf to X0???
	try {
		opencl.error = opencl.getQueue().enqueueCopyBuffer(d_xf, d_x0, 0, 0, realSize * x0elements);
	}
	catch (cl::Error &er) {
		printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
		throw er;
	}    
}


//set new Pars. Cannot update nPts
void CLODE::setPars( std::vector<double>  newPars) {
    
    if (newPars.size()==(size_t)nPts*nPar) {
        pars=newPars;
    
        //sync to device
        try {
			if (clSinglePrecision) { //downcast to float if desired
				std::vector<float> parsF(pars.begin(),pars.end());
				opencl.error = copy(opencl.getQueue(), parsF.begin(), parsF.end(), d_pars);
			}
			else {
				opencl.error = copy(opencl.getQueue(), pars.begin(), pars.end(), d_pars);
			}
        }
        catch (cl::Error &er) {
            printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
            throw er;
        }        
    }
    else {
        printf("Invalid parameter vector: Expected %d*%d elements, recieved %lu\n", nPts, nPar, newPars.size());
        printf("...Parameters were not updated!\n");
        //~ throw std::invalid_argument("Parameter vector has incorrect size.");
    }
}

   
   
void CLODE::setSolverParams(SolverParams<double> newSp){
    
    sp=newSp;
    try {
		if (clSinglePrecision) { //downcast to float if desired
			SolverParams<float> spF=solverParamsToFloat(sp);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_sp, CL_TRUE, 0, sizeof(spF), &spF);
		}
		else {
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_sp, CL_TRUE, 0, sizeof(sp), &sp);
		}
    }
    catch (cl::Error &er) {
        printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
        throw er;
    }
}


//TODO: define an assignment/type cast operator in the struct?
SolverParams<float> CLODE::solverParamsToFloat(SolverParams<double> sp) {
	SolverParams<float> spF;
	spF.dt=sp.dt;
	spF.dtmax=sp.dtmax;
	spF.abstol=sp.abstol;
	spF.reltol=sp.reltol;
	spF.max_steps=sp.max_steps;
	spF.max_store=sp.max_store;
	spF.nout=sp.nout;
	
	return spF;
}
    
    

//populate the RNGstate vector on the device. nPts must be set
void CLODE::seedRNG() {
	//TODO: what is correct method??? here, using MT to get (nRNGstate x nPts) 64bit words
	
	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_int_distribution<unsigned long long> dis;
  
	for (int i=0; i<nRNGstate*nPts; ++i) {
		//~ uint64_t seed = (uint64_t(i) << 32) | i;
		RNGstate[i]=dis(gen);
	}
	
	try {
		opencl.error = copy(opencl.getQueue(), RNGstate.begin(), RNGstate.end(), d_RNGstate);
		
	}
	catch (cl::Error &er) {
		printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
		throw er;
	}
	
}


//populate the RNGstate vector on the device. nPts must be set
void CLODE::seedRNG(int mySeed) {
	
	for (int i=0; i<nRNGstate*nPts; ++i) {
		RNGstate[i]=mySeed+i;
	}
	
	try {
		opencl.error = copy(opencl.getQueue(), RNGstate.begin(), RNGstate.end(), d_RNGstate);
		
	}
	catch (cl::Error &er) {
		printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str() );
		throw er;
	}
	
}


//Simulation routine
void CLODE::transient() {
    
    if (clInitialized) {
		try {
			//kernel args
			cl_transient.setArg(0, d_tspan);
			cl_transient.setArg(1, d_x0); 
			cl_transient.setArg(2, d_pars); 
			cl_transient.setArg(3, d_sp);
			cl_transient.setArg(4, d_xf); 
			cl_transient.setArg(5, d_auxf);
			cl_transient.setArg(6, d_RNGstate); 
			
			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_transient, cl::NullRange, cl::NDRange(nPts), cl::NullRange);
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




std::vector<double> CLODE::getX0() {

		if (clSinglePrecision) { //cast back to double
			std::vector<float> x0F(x0elements);
			opencl.error = copy(opencl.getQueue(), d_x0, x0F.begin(), x0F.end());
			x0.assign(x0F.begin(),x0F.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_x0, x0.begin(), x0.end());
		}
		
		return x0;
}



std::vector<double> CLODE::getXf() {

		if (clSinglePrecision) { //cast back to double
			std::vector<float> xfF(x0elements);
			opencl.error = copy(opencl.getQueue(), d_xf, xfF.begin(), xfF.end());
			xf.assign(xfF.begin(),xfF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_xf, xf.begin(), xf.end());
		}
		
		return xf;
}


std::vector<double> CLODE::getAuxf() { //cast back to double

		if (clSinglePrecision) { 
			std::vector<float> auxfF(auxfelements);
			opencl.error = copy(opencl.getQueue(), d_auxf, auxfF.begin(), auxfF.end());
			auxf.assign(auxfF.begin(),auxfF.end());
		}
		else {
			opencl.error = copy(opencl.getQueue(), d_auxf, auxf.begin(), auxf.end());
		}
		
		return auxf;
}


void CLODE::printStatus() {
	
    opencl.print();
    printf("------------------\n");
    printf("   %s\n", clRHSfilename.c_str());
    printf("   nVar=%d\n", nVar);
    printf("   nPar=%d\n", nPar);
    printf("   nAux=%d\n", nAux);
    printf("   nWiener=%d\n", nWiener);    
	printf("Using %s precision.\n",(clSinglePrecision?"single":"double"));
    printf("Using stepper: %s \n", stepperName.c_str());
}
