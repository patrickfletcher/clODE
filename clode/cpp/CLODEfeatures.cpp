#include "CLODEfeatures.hpp"

#include <algorithm> //std::max
#include <cmath>
#include <stdexcept>

#include "spdlog/spdlog.h"

CLODEfeatures::CLODEfeatures(ProblemInfo prob, std::string stepper, std::string observer, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot)
	: CLODE(prob, stepper, clSinglePrecision, opencl, clodeRoot), observer(observer)
{
	// default fVarIx and eVarIx to allow first query of observer define map (exposes availableObservers, fNames)
	op.fVarIx=0;
	op.eVarIx=0;
	updateObserverDefineMap();

	clprogramstring += read_file(clodeRoot + "initializeObserver.cl");
	clprogramstring += read_file(clodeRoot + "features.cl");
	spdlog::debug("constructor clODEfeatures\n");
}

CLODEfeatures::CLODEfeatures(ProblemInfo prob, std::string stepper, std::string observer, bool clSinglePrecision,  unsigned int platformID, unsigned int deviceID, const std::string clodeRoot)
	: CLODE(prob, stepper, clSinglePrecision, platformID, deviceID, clodeRoot), observer(observer)
{
	// default fVarIx and eVarIx to allow first query of observer define map
	op.fVarIx=0;
	op.eVarIx=0;
	updateObserverDefineMap();

	clprogramstring += read_file(clodeRoot + "initializeObserver.cl");
	clprogramstring += read_file(clodeRoot + "features.cl");
	spdlog::debug("constructor clODEfeatures\n");
}

CLODEfeatures::~CLODEfeatures() {}

// build program and create kernel objects - requires host variables to be set (specifically observerBuildOpts)
void CLODEfeatures::buildCL()
{
	observerBuildOpts=" -D" + observerDefineMap.at(observer).define;
	buildProgram(observerBuildOpts);

	//set up the kernels
	try
	{
		cl_transient = cl::Kernel(opencl.getProgram(), "transient", &opencl.error);
		cl_initializeObserver = cl::Kernel(opencl.getProgram(), "initializeObserver", &opencl.error);
		cl_features = cl::Kernel(opencl.getProgram(), "features", &opencl.error);

		// size_t preferred_multiple;
		// cl::Device dev;
		// opencl.getProgram().getInfo(CL_PROGRAM_DEVICES,&dev);
		// cl_features.getWorkGroupInfo(dev,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,&preferred_multiple);
		// spdlog::info("Preferred work size multiple (features): {}\n",preferred_multiple);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODEfeatures::initializeFeaturesKernel:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("initialize features kernel\n");

	clInitialized = false;
	spdlog::info("Using observer: {}\n",observer.c_str());
}


const std::string CLODEfeatures::getProgramString()
{
	observerBuildOpts=" -D" + observerDefineMap.at(observer).define;
	setCLbuildOpts(observerBuildOpts);
	return buildOptions+clprogramstring+ODEsystemsource; 
}

void CLODEfeatures::initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp) {
    throw std::runtime_error("CLODEfeatures requires the newOp parameter!");
}

//initialize everything
void CLODEfeatures::initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp, ObserverParams<cl_double> newOp)
{
	clInitialized = false;

	//at the time of initialize, make sure observerDataSize and nFeatures are up to date (for d_F, d_odata)
	updateObserverDefineMap();

	setTspan(newTspan);
	setProblemData(newX0, newPars); //will set nPts
	resizeFeaturesVariables(); //set up d_F and d_odata too, which depend on nPts
	setSolverParams(newSp);
	setObserverParams(newOp);

	doObserverInitialization = true;
	clInitialized = true;
	spdlog::debug("initialize clODEfeatures.\n");
}

void CLODEfeatures::setObserver(std::string newObserver)
{
	auto loc = observerDefineMap.find(newObserver); //from steppers.cl
	if ( loc != observerDefineMap.end() )
	{
		observer = newObserver;
		updateObserverDefineMap();
		clInitialized = false;
	}
	else
	{
		spdlog::warn("unknown observer: {}. Observer method unchanged\n",newObserver.c_str());
	}
	spdlog::debug("set observer\n");
}


void CLODEfeatures::setObserverParams(ObserverParams<cl_double> newOp)
{
	try
	{
		if (!clInitialized)
		{
			if (clSinglePrecision)
				d_op = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(ObserverParams<cl_float>), NULL, &opencl.error);
			else
				d_op = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(ObserverParams<cl_double>), NULL, &opencl.error);
		}

		op = newOp; 

		if (clSinglePrecision)
		{ //downcast to float if desired
			ObserverParams<cl_float> opF = observerParamsToFloat(newOp);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_op, CL_TRUE, 0, sizeof(opF), &opF);
		}
		else
		{
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_op, CL_TRUE, 0, sizeof(op), &newOp);
		}

		//if op.fVarIx or op.eVarIx change, observer's fNames may change (doesn't need rebuild though)
		updateObserverDefineMap();
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODEfeatures::setObserverParams:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set observer params\n");
}


// updates the defineMap to reflect changes in problem, precision, observerParams
void CLODEfeatures::updateObserverDefineMap()
{
	getObserverDefineMap(prob, op.fVarIx, op.eVarIx, observerDefineMap, availableObserverNames);
	observerBuildOpts=" -D" + observerDefineMap.at(observer).define;
	if (clSinglePrecision)
		observerDataSize=observerDefineMap.at(observer).observerDataSizeFloat;
	else
		observerDataSize=observerDefineMap.at(observer).observerDataSizeDouble;
	spdlog::debug("observerDataSize = {}\n",observerDataSize);
	observerDataSize = observerDataSize + observerDataSize % realSize; //align to a multiple of realsize. is this necessary?

	nFeatures=(int)observerDefineMap.at(observer).featureNames.size();
	featureNames=observerDefineMap.at(observer).featureNames;
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
		spdlog::info("nPts is too large, requested memory size exceeds selected device's limit. Maximum nPts appears to be {} \n", maxNpts);
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
			spdlog::error("CLODEfeatures::resizeFeaturesVariables:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("resize F, d_F, d_odata with: nPts={}, nF={}\n",nPts,nFeatures);
	}
	
}

//Simulation routines

//
void CLODEfeatures::initializeObserver()
{

	if (clInitialized)
	{
		// spdlog::info("do init={}\n",doObserverInitialization?"true":"false");
		//resize output variables - will only occur if nPts has changed
		resizeFeaturesVariables();

		try
		{
			//kernel arguments
			int ix = 0;
			cl_initializeObserver.setArg(ix++, d_tspan);
			cl_initializeObserver.setArg(ix++, d_x0);
			cl_initializeObserver.setArg(ix++, d_pars);
			cl_initializeObserver.setArg(ix++, d_sp);
			cl_initializeObserver.setArg(ix++, d_RNGstate);
			cl_initializeObserver.setArg(ix++, d_dt);
			cl_initializeObserver.setArg(ix++, d_odata);
			cl_initializeObserver.setArg(ix++, d_op);

			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_initializeObserver, cl::NullRange, cl::NDRange(nPts));
			// spdlog::info("Enqueue error code: {}\n",CLErrorString(opencl.error).c_str());
			opencl.error = opencl.getQueue().finish();
			// spdlog::info("Finish Queue error code: {}\n",CLErrorString(opencl.error).c_str());
			doObserverInitialization = false;
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODEfeatures::initializeObserver:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("run initializeObserver\n");
	}
	else
	{
		spdlog::info("CLODE has not been initialized\n");
	}
}

//overload to allow manual re-initialization of observer data at any time.
void CLODEfeatures::features(bool newDoObserverInitFlag)
{
	doObserverInitialization = newDoObserverInitFlag;

	features();
}

void CLODEfeatures::features()
{

	if (clInitialized)
	{
		// spdlog::info("do init={}\n",doObserverInitialization?"true":"false");
		//resize output variables - will only occur if nPts has changed
		resizeFeaturesVariables();

		if (doObserverInitialization)
			initializeObserver();

		try
		{
			//kernel arguments
			int ix = 0;
			cl_features.setArg(ix++, d_tspan);
			cl_features.setArg(ix++, d_x0);
			cl_features.setArg(ix++, d_pars);
			cl_features.setArg(ix++, d_sp);
			cl_features.setArg(ix++, d_xf);
			cl_features.setArg(ix++, d_RNGstate);
			cl_features.setArg(ix++, d_dt);
			cl_features.setArg(ix++, d_odata);
			cl_features.setArg(ix++, d_op);
			cl_features.setArg(ix++, d_F);

			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_features, cl::NullRange, cl::NDRange(nPts));
			// spdlog::info("Enqueue error code: {}\n",CLErrorString(opencl.error).c_str());
			opencl.error = opencl.getQueue().finish();
			// spdlog::info("Finish Queue error code: {}\n",CLErrorString(opencl.error).c_str());
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODEfeatures::features:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("run features\n");
	}
	else
	{
		spdlog::info("CLODE has not been initialized\n");
	}
}

const std::vector<cl_double> CLODEfeatures::getF()
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