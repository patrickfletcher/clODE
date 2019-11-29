#include "CLODEfeatures.hpp"

// #define dbg_printf printf
#define dbg_printf
#ifdef MATLAB_MEX_FILE
#include "mex.h"
#define printf mexPrintf
#endif

#include <algorithm> //std::max
#include <cmath>
#include <stdexcept>
#include <stdio.h>

CLODEfeatures::CLODEfeatures(ProblemInfo prob, std::string stepper, std::string observer, bool clSinglePrecision, OpenCLResource opencl)
	: CLODE(prob, stepper, clSinglePrecision, opencl)
{

	//printf("\nCLODE has been specialized: CLODEfeatures\n");

	getObserverDefineMap(prob, clSinglePrecision, 2, 0, observerDefineMap, availableObserverNames);
	setObserver(observer);

	clprogramstring += read_file(clodeRoot + "initializeObserver.cl");
	clprogramstring += read_file(clodeRoot + "features.cl");
	dbg_printf("constructor clODEfeatures\n");
}

CLODEfeatures::~CLODEfeatures() {}

//initialize everything
void CLODEfeatures::initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp, ObserverParams<cl_double> newOp)
{
	clInitialized = false;
	// observerBuildOpts = getObserverBuildOpts();
	getObserverDefineMap(prob, clSinglePrecision, newOp.fVarIx, newOp.eVarIx, observerDefineMap, availableObserverNames);
	setObserver(observer);

	//(re)build the program
	buildProgram(observerBuildOpts);

	//Base class initialize kernel:
	initializeTransientKernel();

	//initialize the trajectory kernel
	initializeFeaturesKernel();
	setSolverParams(newSp);
	setObserverParams(newOp);
	setTspan(newTspan);
	setProblemData(newX0, newPars); //will set nPts
	resizeFeaturesVariables(); //set up d_F and d_odata too, which depend on nPts

	printf("Using observer: %s\n",observer);

	clInitialized = true;
	dbg_printf("initialize clODEfeatures.\n");
}

void CLODEfeatures::setObserver(std::string newObserver)
{
	if (newObserver!=observer)
	{
		auto loc = observerDefineMap.find(newObserver); //from steppers.cl
		if ( loc != observerDefineMap.end() )
		{
			observer = newObserver;
			observerBuildOpts=" -D" + observerDefineMap.at(observer).define;
			observerDataSize=observerDefineMap.at(observer).observerDataSize;
			dbg_printf("observerDataSize = %d\n",observerDataSize);
			observerDataSize = observerDataSize + observerDataSize % realSize;
			nFeatures=(int)observerDefineMap.at(observer).featureNames.size();
			featureNames=observerDefineMap.at(observer).featureNames;
			clInitialized = false;
		}
		else
		{
			printf("Warning: unknown observer: %s. Observer method unchanged\n",newObserver.c_str());
		}
		dbg_printf("set observer\n");
	}
	
}

void CLODEfeatures::initializeFeaturesKernel()
{

	try
	{
		if (clSinglePrecision)
			d_op = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(ObserverParams<cl_float>), NULL, &opencl.error);
		else
			d_op = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(ObserverParams<cl_double>), NULL, &opencl.error);

		cl_initializeObserver = cl::Kernel(opencl.getProgram(), "initializeObserver", &opencl.error);
		cl_features = cl::Kernel(opencl.getProgram(), "features", &opencl.error);

		// size_t preferred_multiple;
		// cl::Device dev;
		// opencl.getProgram().getInfo(CL_PROGRAM_DEVICES,&dev);
		// cl_features.getWorkGroupInfo(dev,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,&preferred_multiple);
		// printf("Preferred work size multiple (features): %d\n",preferred_multiple);
	}
	catch (cl::Error &er)
	{
		printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	dbg_printf("initialize features kernel\n");
}

void CLODEfeatures::setObserverParams(ObserverParams<cl_double> newOp)
{

	op = newOp;
	try
	{
		if (clSinglePrecision)
		{ //downcast to float if desired
			ObserverParams<cl_float> opF = observerParamsToFloat(op);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_op, CL_TRUE, 0, sizeof(opF), &opF);
		}
		else
		{
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_op, CL_TRUE, 0, sizeof(op), &op);
		}
	}
	catch (cl::Error &er)
	{
		printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	dbg_printf("set observer params\n");
}

//TODO: define an assignment/type cast operator in the struct?
ObserverParams<cl_float> CLODEfeatures::observerParamsToFloat(ObserverParams<cl_double> op)
{
	ObserverParams<cl_float> opF;

	opF.eVarIx = op.eVarIx;
	opF.fVarIx = op.fVarIx;
	opF.maxEventCount = op.maxEventCount;
	opF.minXamp = op.minXamp;
	opF.nHoodRadius = op.nHoodRadius;
	opF.xUpThresh = op.xUpThresh;
	opF.xDownThresh = op.xDownThresh;
	opF.dxUpThresh = op.dxUpThresh;
	opF.dxDownThresh = op.dxDownThresh;
	opF.eps_dx = op.eps_dx;

	return opF;
}

void CLODEfeatures::resizeFeaturesVariables()
{
	size_t currentFelements = nFeatures * nPts;
	size_t largestAlloc = std::max(nFeatures * realSize, observerDataSize) * nPts;

	if (largestAlloc > opencl.getMaxMemAllocSize())
	{
		int maxNpts = floor(opencl.getMaxMemAllocSize() / (cl_ulong)std::max(nFeatures * realSize, observerDataSize));
		printf("nPts is too large, requested memory size exceeds selected device's limit. Maximum nPts appears to be %d \n", maxNpts);
		throw std::invalid_argument("nPts is too large");
	}

	//resize device variables if nPts changed, or if not yet initialized
	if (!clInitialized || Felements != currentFelements)
	{

		Felements = currentFelements;
		F.resize(currentFelements);

		//resize device variables
		try
		{
			d_odata = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, observerDataSize * nPts, NULL, &opencl.error);
			d_F = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * currentFelements, NULL, &opencl.error);
		}
		catch (cl::Error &er)
		{
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		dbg_printf("resize F, d_F, d_odata with: nPts=%d, nF=%d\n",nPts,nFeatures);
	}
	
}

//Simulation routines

//overload to allow manual re-initialization of observer data at any time.
void CLODEfeatures::features(bool newDoObserverInitFlag)
{
	doObserverInitialization = newDoObserverInitFlag;

	features();
}

//
void CLODEfeatures::initializeObserver()
{

	if (clInitialized)
	{
		// printf("do init=%s\n",doObserverInitialization?"true":"false");
		//resize output variables - will only occur if nPts has changed
		resizeFeaturesVariables();

		try
		{
			//kernel arguments
			cl_initializeObserver.setArg(0, d_tspan);
			cl_initializeObserver.setArg(1, d_x0);
			cl_initializeObserver.setArg(2, d_pars);
			cl_initializeObserver.setArg(3, d_sp);
			cl_initializeObserver.setArg(4, d_RNGstate);
			cl_initializeObserver.setArg(5, d_odata);
			cl_initializeObserver.setArg(6, d_op);

			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_initializeObserver, cl::NullRange, cl::NDRange(nPts));
			// printf("Enqueue error code: %s\n",CLErrorString(opencl.error).c_str());
			opencl.error = opencl.getQueue().finish();
			// printf("Finish Queue error code: %s\n",CLErrorString(opencl.error).c_str());
			doObserverInitialization = 0;
		}
		catch (cl::Error &er)
		{
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		dbg_printf("run initializeObserver\n");
	}
	else
	{
		printf("CLODE has not been initialized\n");
	}
}

void CLODEfeatures::features()
{

	if (clInitialized)
	{
		// printf("do init=%s\n",doObserverInitialization?"true":"false");
		//resize output variables - will only occur if nPts has changed
		resizeFeaturesVariables();

		if (doObserverInitialization)
			initializeObserver();

		try
		{
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

			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_features, cl::NullRange, cl::NDRange(nPts));
			// printf("Enqueue error code: %s\n",CLErrorString(opencl.error).c_str());
			opencl.error = opencl.getQueue().finish();
			// printf("Finish Queue error code: %s\n",CLErrorString(opencl.error).c_str());
		}
		catch (cl::Error &er)
		{
			printf("ERROR: %s(%s)\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		dbg_printf("run features\n");
	}
	else
	{
		printf("CLODE has not been initialized\n");
	}
}

std::vector<double> CLODEfeatures::getF()
{

	if (clSinglePrecision)
	{ //cast back to double
		std::vector<cl_float> FF(Felements);
		opencl.error = copy(opencl.getQueue(), d_F, FF.begin(), FF.end());
		F.assign(FF.begin(), FF.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_F, F.begin(), F.end());
	}

	return F;
}
