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
	getStepperDefineMap(stepperDefineMap, availableSteppers); // from steppers.cl
	setProblemInfo(prob);
	setStepper(stepper);
	setPrecision(clSinglePrecision);
	setClodeRoot(clodeRoot);
	setOpenCL(opencl);

	clprogramstring = read_file(clodeRoot + "transient.cl");
	spdlog::debug("constructor clODE");
}

CLODE::CLODE(ProblemInfo prob, std::string stepper, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot)
{
	getStepperDefineMap(stepperDefineMap, availableSteppers); // from steppers.cl
	setProblemInfo(prob);
	setStepper(stepper);
	setPrecision(clSinglePrecision);
	setClodeRoot(clodeRoot);
	setOpenCL(platformID, deviceID);

	clprogramstring = read_file(clodeRoot + "transient.cl");
	spdlog::debug("constructor clODE");
}

CLODE::~CLODE()
{
}

void CLODE::setProblemInfo(ProblemInfo newProb)
{ // TODO: not equality check for ProblemInfo struct, error checking: at least one variable!
	prob = newProb;
	clRHSfilename = newProb.clRHSfilename;
	ODEsystemsource = read_file(clRHSfilename);
	nVar = newProb.nVar;
	nPar = newProb.nPar;
	nAux = newProb.nAux;
	nWiener = newProb.nWiener;

	// Some OpenCL implementations don't support zero-length arrays in OpenCL programs
	// nPar = nPar>0?nPar:1;
	// nAux = nAux>0?nAux:1;
	// nWiener = nWiener>0?nWiener:1;

	spdlog::debug("set new problem");
}

void CLODE::setStepper(std::string newStepper)
{
	// if (newStepper!=stepper)
	// {
	auto loc = stepperDefineMap.find(newStepper); // from steppers.cl
	if (loc != stepperDefineMap.end())
	{
		stepper = newStepper;
	}
	else
	{
		spdlog::warn("Unknown stepper: {}. Stepper method unchanged", newStepper.c_str());
	}
	spdlog::debug("set stepper");
	// }
}

void CLODE::setPrecision(bool newPrecision)
{
	// if (newPrecision != clSinglePrecision)
	// {
	clSinglePrecision = newPrecision;
	realSize = newPrecision ? sizeof(cl_float) : sizeof(cl_double);
	spdlog::debug("set precision");
	// }
}

void CLODE::setOpenCL(OpenCLResource newOpencl)
{ // TODO: not equality check for OpenCLResource class
	//~ if (newOpencl!=opencl) {
	opencl = newOpencl;
	//~ }
	spdlog::debug("set OpenCL");
}

void CLODE::setOpenCL(unsigned int platformID, unsigned int deviceID)
{ // TODO: not equality check for OpenCLResource class
	//~ if (newOpencl!=opencl) {
	opencl = OpenCLResource(platformID, deviceID);
	//~ }
	spdlog::debug("set OpenCL");
}

void CLODE::setClodeRoot(const std::string newClodeRoot)
{
	clodeRoot = newClodeRoot;
}

void CLODE::setCLbuildOpts(std::string extraBuildOpts)
{

	if (!clSinglePrecision && !opencl.getDoubleSupport())
	{ // TODO: make this an error?
		clSinglePrecision = true;
		spdlog::warn("device selected does not support double precision. Using single precision");
	}

	buildOptions = "";

	// specify precision
	if (clSinglePrecision)
		buildOptions += " -DCLODE_SINGLE_PRECISION";
	else
		buildOptions += " -DCLODE_DOUBLE_PRECISION";

	// specify stepper
	buildOptions += " -D" + stepperDefineMap.at(stepper);

	// specify problem dimensions
	buildOptions += " -DN_PAR=" + std::to_string((long long)nPar); // for older c++ compilers the to_string(int) overload of the STL isn't present
	buildOptions += " -DN_VAR=" + std::to_string((long long)nVar);
	buildOptions += " -DN_AUX=" + std::to_string((long long)nAux);
	buildOptions += " -DN_WIENER=" + std::to_string((long long)nWiener);

	// include folder for CLODE
	buildOptions += " -I" + clodeRoot;

	buildOptions += extraBuildOpts;

	spdlog::debug("OpenCL build options {}", buildOptions);
}

// build creates build option defined constants based on selected options, adds the ODEsystem source to clprogramstring then builds for selected OpenCL resource
void CLODE::buildProgram(std::string extraBuildOpts)
{
	setCLbuildOpts(extraBuildOpts);

	opencl.buildProgramFromString(clprogramstring + ODEsystemsource, buildOptions);

	spdlog::trace(clprogramstring + ODEsystemsource);
	spdlog::debug("CLODE buildProgram finished");
}

// build program and create kernel objects. requires host variables to be set
void CLODE::buildCL()
{
	spdlog::info("Running CLODE buildCL");
	buildProgram();

	// set up the kernel
	try
	{
		cl_transient = cl::Kernel(opencl.getProgram(), "transient", &opencl.error);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::buildCL(): create kernels{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("Created kernel");
}

// resize all the nPts dependent variables, only if nPts changed
void CLODE::setNpts(cl_int newNpts)
{
	size_t largestAlloc = std::max(nVar, std::max(nPar, nAux)) * nPts * realSize;
	// spdlog::info("Computed largestAlloc: {}", largestAlloc);

	if (largestAlloc > opencl.getMaxMemAllocSize())
	{
		throw std::invalid_argument("nPts*nVar, nPts*nPar, or nPts*nAux is too large");
	}

	if (newNpts != nPts)
	{
		nPts = newNpts;

		x0elements = nVar * nPts;
		parselements = nPar * nPts;

		RNGelements = nRNGstate * nPts;
		// cl::Buffer doesn't like zero-sized arrays:
		// RNGelements = nWiener > 0 ? nWiener * nRNGstate * nPts : 1;

		// resize host variables
		x0.resize(x0elements);
		pars.resize(parselements);

		RNGstate.resize(RNGelements);

		dt.resize(nPts);
		std::fill(dt.begin(), dt.end(), sp.dt);

		xf.resize(x0elements);
		// new device variables
		try
		{
			// must be populated with valid data prior to simulation
			d_x0 = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * x0elements, NULL, &opencl.error);
			d_pars = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, realSize * parselements, NULL, &opencl.error);

			d_RNGstate = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, sizeof(cl_ulong) * RNGelements, NULL, &opencl.error);

			// resize and fill dt buffer
			d_dt = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * nPts, NULL, &opencl.error);
			if (clSinglePrecision) {
				std::vector<cl_float> dtF(dt.begin(), dt.end());
				opencl.error = copy(opencl.getQueue(), dtF.begin(), dtF.end(), d_dt);
			}
			else
				opencl.error = copy(opencl.getQueue(), dt.begin(), dt.end(), d_dt);

			// output of simulation
			d_xf = cl::Buffer(opencl.getContext(), CL_MEM_READ_WRITE, realSize * x0elements, NULL, &opencl.error);
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODE::setNpts:{}({})", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}

		// seed RNG must occur after device variable d_RNGstate is resized
		seedRNG();
		spdlog::debug("set nPts={}", nPts);
	}
}

// initialize new set of trajectories (nPts may change)
void CLODE::setProblemData(std::vector<cl_double> newX0, std::vector<cl_double> newPars)
{ // check if newX0 and newPars are valid, and update nPts if needed:
	if (newX0.size() % nVar != 0)
	{
		spdlog::info("Invalid initial condition vector: not a multiple of nVar={}", nVar);
		spdlog::info("...Initial conditions were not updated!");
		return;
	}

	if (newPars.size() % nPar != 0)
	{
		spdlog::info("Invalid parameter vector: not a multiple of nPar={}", nPar);
		spdlog::info("...Parameters were not updated!");
		return;
	}

	// now check if newX0 and newPars represent same number of sets
	cl_int nPtsX0 = newX0.size() / nVar;
	cl_int nPtsPars = newPars.size() / nPar;
	// spdlog::info("Computed nPts: {} {}", nPtsX0, nPtsPars);
	if (nPtsX0 != nPtsPars)
	{
		spdlog::info("Initial contition and parameter vector dimensions don't match");
		spdlog::info("...Expected {} sets of each, recieved {} for x0 and {} for pars", nPts, nPtsX0, nPtsPars);
		spdlog::info("...Problem data was not updated!");
		return;
	}

	// set nPts
	setNpts(nPtsX0);

	// set things that depend on nPts
	setX0(newX0);
	setPars(newPars);
	spdlog::debug("set problem data");
}

// set new x0. Cannot update nPts
void CLODE::setX0(std::vector<cl_double> newX0)
{
	if (newX0.size() == (size_t)nPts * nVar)
	{
		x0 = newX0;

		// sync to device
		try
		{
			if (clSinglePrecision)
			{ // downcast to float if desired
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
			spdlog::error("CLODE::setX0:{}({})", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("set X0");
	}
	else
	{
		// spdlog::info("Invalid initial condition vector: Expected {}*{} elements, recieved {}}", nPts, nVar, newX0.size());
		spdlog::info("...Initial conditions were not updated!");
		//~ throw std::invalid_argument("Initial Condition vector has incorrect size.");
	}
}

// set new Pars. Cannot update nPts
void CLODE::setPars(std::vector<cl_double> newPars)
{
	if (newPars.size() == (size_t)nPts * nPar)
	{
		pars = newPars;

		// sync to device
		try
		{
			if (clSinglePrecision)
			{ // downcast to float if desired
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
			spdlog::error("CLODE::setPars:{}({})", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("set P");
	}
	else
	{
		spdlog::info("Invalid parameter vector: Expected {}*{} elements, recieved {}", nPts, nPar, newPars.size());
		spdlog::info("...Parameters were not updated!");
		//~ throw std::invalid_argument("Parameter vector has incorrect size.");
	}
}

void CLODE::setTspan(std::vector<cl_double> newTspan)
{
	try
	{
		tspan = newTspan;
		d_tspan = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, realSize * 2, NULL, &opencl.error);

		if (clSinglePrecision)
		{ // downcast to float if desired
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
		spdlog::error("CLODE::setTspan:{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set tspan");
}

void CLODE::setSolverParams(SolverParams<cl_double> newSp)
{ // TODO: equality operator for SolverParams struct
	try
	{
		sp = newSp;
		// std::fill(dt.begin(), dt.end(), sp.dt);

		if (clSinglePrecision)
		{ // downcast to float if desired
			SolverParams<cl_float> spF = solverParamsToFloat(sp);
			d_sp = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(SolverParams<cl_float>), NULL, &opencl.error);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_sp, CL_TRUE, 0, sizeof(spF), &spF);
			// std::vector<cl_float> dtF(dt.begin(), dt.end());
			// opencl.error = copy(opencl.getQueue(), dtF.begin(), dtF.end(), d_dt);
			// opencl.error = opencl.getQueue().enqueueFillBuffer(d_dt, spF.dt, 0, sizeof(spF.dt));
		}
		else
		{
			d_sp = cl::Buffer(opencl.getContext(), CL_MEM_READ_ONLY, sizeof(SolverParams<cl_double>), NULL, &opencl.error);
			opencl.error = opencl.getQueue().enqueueWriteBuffer(d_sp, CL_TRUE, 0, sizeof(sp), &sp);
			// opencl.error = copy(opencl.getQueue(), dt.begin(), dt.end(), d_dt);
			// opencl.error = opencl.getQueue().enqueueFillBuffer(d_dt, sp.dt, 0, sizeof(sp.dt));
		}
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::setSolverParams:{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set SolverParams");
}

// TODO: define an assignment/type cast operator in the struct?
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

// populate the RNGstate vector on the device. nPts must be set
void CLODE::seedRNG()
{
	// TODO: what is correct method??? here, using MT to get (nRNGstate x nPts) 64bit words

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
		spdlog::error("CLODE::seedRNG:{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set random RNG seed");
}

// populate the RNGstate vector on the device. nPts must be set
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
		spdlog::error("CLODE::seedRNG(int mySeed):{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("set fixed RNG seed");
}

// Simulation routine
void CLODE::transient()
{
	try
	{
		// kernel args
		int ix = 0;
		cl_transient.setArg(ix++, d_tspan);
		cl_transient.setArg(ix++, d_x0);
		cl_transient.setArg(ix++, d_pars);
		cl_transient.setArg(ix++, d_sp);
		cl_transient.setArg(ix++, d_xf);
		cl_transient.setArg(ix++, d_RNGstate);
		cl_transient.setArg(ix++, d_dt);

		// execute the kernel
		opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_transient, cl::NullRange, cl::NDRange(nPts));
		opencl.getQueue().finish();
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::transient:{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::info("run transient");
}

void CLODE::shiftTspan()
{
	std::vector<cl_double> newTspan({tspan[1], tspan[1] + (tspan[1] - tspan[0])});
	setTspan(newTspan);
	spdlog::debug("shift tspan");
}

void CLODE::shiftX0()
{
	// device to device transfer of Xf to X0
	try
	{
		opencl.error = opencl.getQueue().enqueueCopyBuffer(d_xf, d_x0, 0, 0, realSize * x0elements);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODE::shiftX0:{}({})", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	spdlog::debug("shift X0");
}

const std::vector<cl_double> CLODE::getX0()
{

	if (clSinglePrecision)
	{ // cast back to double
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

const std::vector<cl_double> CLODE::getXf()
{

	if (clSinglePrecision)
	{ // cast back to double
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

void CLODE::printStatus()
{
	// opencl.print();
	spdlog::info("------------------");
	spdlog::info("   {}", clRHSfilename.c_str());
	spdlog::info("   nVar={}", nVar);
	spdlog::info("   nPar={}", nPar);
	spdlog::info("   nAux={}", nAux);
	spdlog::info("   nWiener={}", nWiener);
	spdlog::info("Using {} precision.", (clSinglePrecision ? "single" : "double"));
	spdlog::info("Using stepper: {} ", stepper.c_str());
}
