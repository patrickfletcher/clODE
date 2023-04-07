#include "CLODE.hpp"
#include "clODE_struct_defs.cl"
#include "OpenCLResource.hpp"
#include "steppers.cl"

#include <algorithm> //std::max
#include <cmath>
#include <random>
#include <stdexcept>

#include "spdlog/spdlog.h"

CLODE::CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot)
{
	getStepperDefineMap(stepperDefineMap, availableSteppers); //from steppers.cl
	setNewProblem(prob);
	setStepper(stepper);
	setPrecision(clSinglePrecision);
    setClodeRoot(clodeRoot);
	setOpenCL(opencl);

	clprogramstring = read_file(clodeRoot + "transient.cl");
	spdlog::debug("constructor clODE\n");
}

CLODE::CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot)
{
	getStepperDefineMap(stepperDefineMap, availableSteppers); //from steppers.cl
	setNewProblem(prob);
	setStepper(stepper);
	setPrecision(clSinglePrecision);
    setClodeRoot(clodeRoot);
	setOpenCL(platformID, deviceID);

	clprogramstring = read_file(clodeRoot + "transient.cl");
	spdlog::debug("constructor clODE\n");
}

CLODE::~CLODE()
{
}

void CLODE::setNewProblem(ProblemInfo newProb)
{ //TODO: not equality check for ProblemInfo struct, error checking: at least one variable!
	prob=newProb;
	clRHSfilename = newProb.clRHSfilename;
	ODEsystemsource = read_file(clRHSfilename);
	nVar = newProb.nVar;
	nPar = newProb.nPar>0?newProb.nPar:1; //support zero params
	nAux = newProb.nAux>0?newProb.nAux:1; //support zero aux
	nWiener = newProb.nWiener>0?newProb.nWiener:1; //support zero wiener

	clInitialized = false;
	spdlog::debug("set new problem\n");
}

void CLODE::setStepper(std::string newStepper)
{
	// if (newStepper!=stepper)
	// {
	auto loc = stepperDefineMap.find(newStepper); //from steppers.cl
	if ( loc != stepperDefineMap.end() )
	{
		stepper = newStepper;
		clInitialized = false;
	}
	else
	{
		spdlog::warn("Unknown stepper: {}. Stepper method unchanged\n",newStepper.c_str());
	}
	spdlog::debug("set stepper\n");	
	// }
}

void CLODE::setPrecision(bool newPrecision)
{
	// if (newPrecision != clSinglePrecision)
	// {
	clSinglePrecision = newPrecision;
	realSize = newPrecision ? sizeof(cl_float) : sizeof(cl_double);
	clInitialized = false;
	spdlog::debug("set precision\n");
	// }
}

void CLODE::setOpenCL(OpenCLResource newOpencl)
{//TODO: not equality check for OpenCLResource class
	//~ if (newOpencl!=opencl) {
	opencl = newOpencl;
	clInitialized = false;
	//~ }
	spdlog::debug("set OpenCL\n");
}

void CLODE::setOpenCL(unsigned int platformID, unsigned int deviceID)
{//TODO: not equality check for OpenCLResource class
	//~ if (newOpencl!=opencl) {
	opencl = OpenCLResource(platformID, deviceID);
	clInitialized = false;
	//~ }
	spdlog::debug("set OpenCL\n");
}

void CLODE::setClodeRoot(const std::string newClodeRoot)
{
    clodeRoot = newClodeRoot;
}


void CLODE::setCLbuildOpts(std::string extraBuildOpts)
{
	
	if (!clSinglePrecision && !opencl.getDoubleSupport())
	{ //TODO: make this an error?
		clSinglePrecision = true;
		spdlog::warn("device selected does not support double precision. Using single precision\n");
	}

	buildOptions = "";

	//specify precision
	if (clSinglePrecision)
		buildOptions += " -DCLODE_SINGLE_PRECISION";
	else
		buildOptions += " -DCLODE_DOUBLE_PRECISION";

	// spdlog::info("{}\n",stepperDefineMap.at(stepper).c_str());

	//specify stepper
	buildOptions += " -D" + stepperDefineMap.at(stepper);
	// buildOptions += getStepperDefine();

	//specify problem dimensions
	buildOptions += " -DN_PAR=" + std::to_string((long long)nPar); //for older c++ compilers the to_string(int) overload of the STL isn't present
	buildOptions += " -DN_VAR=" + std::to_string((long long)nVar);
	buildOptions += " -DN_AUX=" + std::to_string((long long)nAux);
	buildOptions += " -DN_WIENER=" + std::to_string((long long)nWiener);

	//include folder for CLODE
	buildOptions += " -I" + clodeRoot;
    //buildOptions += " -I" + clodeRoot + "../../";

    spdlog::debug("OpenCL build options {}", buildOptions);

	buildOptions += extraBuildOpts;
}


//build creates build option defined constants based on selected options, adds the ODEsystem source to clprogramstring then builds for selected OpenCL resource
void CLODE::buildProgram(std::string extraBuildOpts)
{
	setCLbuildOpts(extraBuildOpts);

	//ODEsystem source is delayed to here, in case we change it
	// ODEsystemsource = read_file(clRHSfilename);
	// clprogramstring += ODEsystemsource;

	// spdlog::info("{}", clprogramstring.c_str());
	// spdlog::info("{}", ODEsystemsource.c_str());
	// spdlog::info("{}", buildOptions.c_str());

	//now build

    spdlog::debug(clprogramstring + ODEsystemsource);
    spdlog::debug(buildOptions);
	opencl.buildProgramFromString(clprogramstring + ODEsystemsource, buildOptions);

	// printStatus();
	spdlog::debug("build clODE\n");
}

// build program and create kernel objects. requires host variables to be set
void CLODE::buildCL()
{
	buildProgram();

	//set up the kernel
	try
	{ 
		cl_transient = cl::Kernel(opencl.getProgram(), "transient", &opencl.error);

		// size_t preferred_multiple;
		// cl::Device dev;
		// opencl.getProgram().getInfo(CL_PROGRAM_DEVICES,&dev);
		// cl_transient.getWorkGroupInfo(dev,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,&preferred_multiple);
		// spdlog::info("Preferred work size multiple (transient): {}\n",preferred_multiple);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::initializeTransientKernel():{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	clInitialized = false;
	spdlog::debug("buildCL\n");
}

//initialize everything: build the program, create the kernels, and set all needed problem data.
void CLODE::initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp)
{
	clInitialized = false;

	setTspan(newTspan);
	setProblemData(newX0, newPars); //will call setNpts
	setSolverParams(newSp);

	clInitialized = true;
	spdlog::debug("initialize clODE\n");
}

//initialize new set of trajectories (nPts may change)
void CLODE::setProblemData(std::vector<cl_double> newX0, std::vector<cl_double> newPars)
{	//check if newX0 and newPars are valid, and update nPts if needed:
	if (newX0.size() % nVar != 0)
	{
		spdlog::info("Invalid initial condition vector: not a multiple of nVar={}\n", nVar);
		spdlog::info("...Initial conditions were not updated!\n");
		return;
	}

	if (newPars.size() % nPar != 0)
	{
		spdlog::info("Invalid parameter vector: not a multiple of nPar={}\n", nPar);
		spdlog::info("...Parameters were not updated!\n");
		return;
	}

	// now check if newX0 and newPars represent same number of sets
	cl_int nPtsX0 = newX0.size() / nVar;
	cl_int nPtsPars = newPars.size() / nPar;
	// spdlog::info("Computed nPts: {} {}\n", nPtsX0, nPtsPars);
	if (nPtsX0 != nPtsPars)
	{
		spdlog::info("Initial contition and parameter vector dimensions don't match");
		spdlog::info("...Expected {} sets of each, recieved {} for x0 and {} for pars\n", nPts, nPtsX0, nPtsPars);
		spdlog::info("...Problem data was not updated!\n");
		return;
	}

	//set nPts
	setNpts(nPtsX0);

	//set things that depend on nPts
	setX0(newX0);
	setPars(newPars);
	spdlog::debug("set problem data\n");
}

//resize all the nPts dependent variables, only if nPts changed
void CLODE::setNpts(cl_int newNpts)
{	//unlikely that any of these should ever exceed memory limits...
	size_t largestAlloc = std::max(nVar, std::max(nPar, nAux)) * nPts * realSize;
	// spdlog::info("Computed largestAlloc: {}\n", largestAlloc);

	if (largestAlloc > opencl.getMaxMemAllocSize())
	{
		throw std::invalid_argument("nPts*nVar, nPts*nPar, or nPts*nAux is too large");
	}

	if (!clInitialized || newNpts != nPts)
	{
		nPts = newNpts;

		x0elements = nVar * nPts;
		parselements = nPar * nPts;
		RNGelements = nRNGstate * nPts;

		//resize host variables
		x0.resize(x0elements);
		pars.resize(parselements);
		xf.resize(x0elements);
		RNGstate.resize(RNGelements);
		dt.resize(nPts);

		//new device variables
		try
		{
			d_x0 = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * x0elements, NULL, &opencl.error);
			d_pars = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, realSize * parselements, NULL, &opencl.error);
			d_xf = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * x0elements, NULL, &opencl.error);
			d_RNGstate = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, sizeof(cl_ulong) * RNGelements, NULL, &opencl.error);
			d_dt = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * nPts, NULL, &opencl.error);
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODE::setNpts:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}

		//seed RNG must occur after device variable d_RNGstate is resized
		seedRNG();
		spdlog::debug("set nPts\n");
	}
}

void CLODE::setTspan(std::vector<cl_double> newTspan)
{
	try
	{
		if (!clInitialized)
			d_tspan = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, realSize * 2, NULL, &opencl.error);

		tspan = newTspan;

		if (clSinglePrecision)
		{ //downcast to float if desired
			std::vector<cl_float> tspanF(tspan.begin(), tspan.end());
			opencl.error = copy(opencl.getQueue(), tspanF.begin(), tspanF.end(), d_tspan);
		}
		else
		{
			opencl.error = copy(opencl.getQueue(), tspan.begin(), tspan.end(), d_tspan);
		}
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::setTspan:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set tspan\n");
}

void CLODE::shiftTspan()
{
	std::vector<cl_double> newTspan({tspan[1], tspan[1] + (tspan[1] - tspan[0])});
	setTspan(newTspan);
	spdlog::debug("shift tspan\n");
}

//set new x0. Cannot update nPts
void CLODE::setX0(std::vector<cl_double> newX0)
{
	if (newX0.size() == (size_t)nPts * nVar)
	{
		x0 = newX0;

		//sync to device
		try
		{
			if (clSinglePrecision)
			{ //downcast to float if desired
				std::vector<cl_float> x0F(x0.begin(), x0.end());
				opencl.error = copy(opencl.getQueue(), x0F.begin(), x0F.end(), d_x0);
			}
			else
			{
				opencl.error = copy(opencl.getQueue(), x0.begin(), x0.end(), d_x0);
			}
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODE::setX0:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("set X0\n");
	}
	else
	{
		// spdlog::info("Invalid initial condition vector: Expected {}*{} elements, recieved {}}\n", nPts, nVar, newX0.size());
		spdlog::info("...Initial conditions were not updated!\n");
		//~ throw std::invalid_argument("Initial Condition vector has incorrect size.");
	}
}

void CLODE::shiftX0()
{
	//device to device transfer of Xf to X0
	try
	{
		opencl.error = opencl.getQueue().enqueueCopyBuffer(d_xf, d_x0, 0, 0, realSize * x0elements);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::shiftX0:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("shift X0\n");
}

//set new Pars. Cannot update nPts
void CLODE::setPars(std::vector<cl_double> newPars)
{
	if (newPars.size() == (size_t)nPts * nPar)
	{
		pars = newPars;

		//sync to device
		try
		{
			if (clSinglePrecision)
			{ //downcast to float if desired
				std::vector<cl_float> parsF(pars.begin(), pars.end());
				opencl.error = copy(opencl.getQueue(), parsF.begin(), parsF.end(), d_pars);
			}
			else
			{
				opencl.error = copy(opencl.getQueue(), pars.begin(), pars.end(), d_pars);
			}
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODE::setPars:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("set P\n");
	}
	else
	{
		spdlog::info("Invalid parameter vector: Expected {}*{} elements, recieved {}\n", nPts, nPar, newPars.size());
		spdlog::info("...Parameters were not updated!\n");
		//~ throw std::invalid_argument("Parameter vector has incorrect size.");
	}
}

void CLODE::setSolverParams(SolverParams<cl_double> newSp)
{//TODO: equality operator for SolverParams struct
	try
	{
		if (!clInitialized)
		{
			if (clSinglePrecision)
				d_sp = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(SolverParams<cl_float>), NULL, &opencl.error);
			else
				d_sp = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(SolverParams<cl_double>), NULL, &opencl.error);	
		}
	
		sp = newSp;
		std::fill(dt.begin(), dt.end(), sp.dt);
		
		if (clSinglePrecision)
		{ //downcast to float if desired
			SolverParams<cl_float> spF = solverParamsToFloat(sp);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_sp, CL_TRUE, 0, sizeof(spF), &spF);
			std::vector<cl_float> dtF(dt.begin(), dt.end());
			opencl.error = copy(opencl.getQueue(), dtF.begin(), dtF.end(), d_dt);
			// opencl.error = opencl.getQueue().enqueueFillBuffer(d_dt, spF.dt, 0, sizeof(spF.dt));
		}
		else
		{
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_sp, CL_TRUE, 0, sizeof(sp), &sp);
			opencl.error = copy(opencl.getQueue(), dt.begin(), dt.end(), d_dt);
			// opencl.error = opencl.getQueue().enqueueFillBuffer(d_dt, sp.dt, 0, sizeof(sp.dt));
		}
		
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::setSolverParams:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set SolverParams\n");
}

//TODO: define an assignment/type cast operator in the struct?
SolverParams<cl_float> CLODE::solverParamsToFloat(SolverParams<cl_double> sp)
{
	SolverParams<cl_float> spF;
	spF.dt = sp.dt;
	spF.dtmax = sp.dtmax;
	spF.abstol = sp.abstol;
	spF.reltol = sp.reltol;
	spF.max_steps = sp.max_steps;
	spF.max_store = sp.max_store;
	spF.nout = sp.nout;

	return spF;
}

//populate the RNGstate vector on the device. nPts must be set
void CLODE::seedRNG()
{
	//TODO: what is correct method??? here, using MT to get (nRNGstate x nPts) 64bit words

	std::random_device rd;
	std::mt19937_64 gen(rd());
	std::uniform_int_distribution<cl_ulong> dis;

	for (int i = 0; i < nRNGstate * nPts; ++i)
	{
		//~ uint64_t seed = (uint64_t(i) << 32) | i;
		RNGstate[i] = dis(gen);
	}

	try
	{
		opencl.error = copy(opencl.getQueue(), RNGstate.begin(), RNGstate.end(), d_RNGstate);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::seedRNG:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set random RNG seed\n");
}

//populate the RNGstate vector on the device. nPts must be set
void CLODE::seedRNG(cl_int mySeed)
{

	for (int i = 0; i < nRNGstate * nPts; ++i)
	{
		RNGstate[i] = mySeed + i;
	}

	try
	{
		opencl.error = copy(opencl.getQueue(), RNGstate.begin(), RNGstate.end(), d_RNGstate);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::seedRNG(int mySeed):{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set fixed RNG seed\n");
}

//Simulation routine
void CLODE::transient()
{

	if (clInitialized)
	{
		try
		{
			//kernel args
			int ix=0;
			cl_transient.setArg(ix++, d_tspan);
			cl_transient.setArg(ix++, d_x0);
			cl_transient.setArg(ix++, d_pars);
			cl_transient.setArg(ix++, d_sp);
			cl_transient.setArg(ix++, d_xf);
			cl_transient.setArg(ix++, d_RNGstate);
			cl_transient.setArg(ix++, d_dt);

			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_transient, cl::NullRange, cl::NDRange(nPts));
			opencl.getQueue().finish();
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODE::transient:{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::info("run transient\n");
	}
	else
	{
		spdlog::info("CLODE has not been initialized\n");
	}
}

std::vector<cl_double> CLODE::getX0()
{

	if (clSinglePrecision)
	{ //cast back to double
		std::vector<cl_float> x0F(x0elements);
		opencl.error = copy(opencl.getQueue(), d_x0, x0F.begin(), x0F.end());
		x0.assign(x0F.begin(), x0F.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_x0, x0.begin(), x0.end());
	}

	return x0;
}

std::vector<cl_double> CLODE::getXf()
{

	if (clSinglePrecision)
	{ //cast back to double
		std::vector<cl_float> xfF(x0elements);
		opencl.error = copy(opencl.getQueue(), d_xf, xfF.begin(), xfF.end());
		xf.assign(xfF.begin(), xfF.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_xf, xf.begin(), xf.end());
	}

	return xf;
}


std::string CLODE::getProgramString() 
{
	setCLbuildOpts();
	return buildOptions+clprogramstring+ODEsystemsource; 
}


void CLODE::printStatus()
{

	opencl.print();
	spdlog::info("------------------\n");
	spdlog::info("   {}\n", clRHSfilename.c_str());
	spdlog::info("   nVar={}\n", nVar);
	spdlog::info("   nPar={}\n", nPar);
	spdlog::info("   nAux={}\n", nAux);
	spdlog::info("   nWiener={}\n", nWiener);
	spdlog::info("Using {} precision.\n", (clSinglePrecision ? "single" : "double"));
	spdlog::info("Using stepper: {} \n", stepper.c_str());
}
