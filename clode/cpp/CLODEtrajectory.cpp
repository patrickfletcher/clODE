#include "CLODEtrajectory.hpp"

#include <algorithm> //std::max
#include <cmath>
#include <stdexcept>

#include "spdlog/spdlog.h"

CLODEtrajectory::CLODEtrajectory(ProblemInfo prob, std::string stepper, bool clSinglePrecision, OpenCLResource opencl, const std::string clodeRoot)
	: CLODE(prob, stepper, clSinglePrecision, opencl, clodeRoot), nStoreMax(0)
{

	//printf("\nCLODE has been specialized: CLODEtrajectory\n");

	clprogramstring += read_file(clodeRoot + "trajectory.cl");
	spdlog::debug("constructor clODEtrajectory\n");
}

CLODEtrajectory::CLODEtrajectory(ProblemInfo prob, std::string stepper, bool clSinglePrecision, unsigned int platformID, unsigned int deviceID, const std::string clodeRoot)
	: CLODE(prob, stepper, clSinglePrecision, platformID, deviceID, clodeRoot), nStoreMax(0)
{

	//printf("\nCLODE has been specialized: CLODEtrajectory\n");

	clprogramstring += read_file(clodeRoot + "trajectory.cl");
	spdlog::debug("constructor clODEtrajectory\n");
}

CLODEtrajectory::~CLODEtrajectory() {}


// build program and create kernel objects. requires host variables to be set
void CLODEtrajectory::buildCL()
{
	buildProgram();

	//set up the kernels
	try
	{
		cl_transient = cl::Kernel(opencl.getProgram(), "transient", &opencl.error);
		cl_trajectory = cl::Kernel(opencl.getProgram(), "trajectory", &opencl.error);

		// size_t preferred_multiple;
		// cl::Device dev;
		// opencl.getProgram().getInfo(CL_PROGRAM_DEVICES,&dev);
		// cl_trajectory.getWorkGroupInfo(dev,CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE,&preferred_multiple);
		// printf("Preferred work size multiple (trajectory): {}\n",preferred_multiple);
	}
	catch (cl::Error &er)
	{
		spdlog::error("CLODEtrajectory::initializeTrajectoryKernel():{}({})\n", er.what(), CLErrorString(er.err()).c_str());
		throw er;
	}
	clInitialized = false;
	spdlog::debug("initialize trajectory kernel\n");
}

//initialize everything
void CLODEtrajectory::initialize(std::vector<cl_double> newTspan, std::vector<cl_double> newX0, std::vector<cl_double> newPars, SolverParams<cl_double> newSp)
{

	clInitialized = false;

	setTspan(newTspan);
	setProblemData(newX0, newPars); //will set nPts
	setSolverParams(newSp);
	//set up output variables, depends on sp.max_store, nPts, nVar, nAux
	resizeTrajectoryVariables(); 

	clInitialized = true;
	spdlog::debug("initialize clODEfeatures\n");
}

void CLODEtrajectory::resizeTrajectoryVariables()
{

	int currentStoreAlloc = sp.max_store; //TODO: always alloc sp.max_store? or try to compute something?
	//catch any changes to tspan, dt, or nout:
	// int currentStoreAlloc=(tspan[1]-tspan[0])/(sp.dt*sp.nout)+1;
	// currentStoreAlloc=std::min(currentStoreAlloc, sp.max_store);
	//~ if (currentStoreAlloc > sp.max_store) {
	//~ currentStoreAlloc=sp.max_store;
	//~ printf("Warning: (tspan[1]-tspan[0])/(dt*nout)+1 > max_store. Reducing storage to max_store={} timepoints. \n",sp.max_store);
	//~ }

	//check largest desired memory chunk against device's maximum allowable variable size
	size_t largestAlloc = std::max(1 /*t*/, std::max(nVar /*x, dx*/, nAux /*aux*/)) * nPts * currentStoreAlloc * realSize;

	if (largestAlloc > opencl.getMaxMemAllocSize())
	{
		int estimatedMaxStoreAlloc = std::floor(opencl.getMaxMemAllocSize() / (std::max(1, std::max(nVar, nAux)) * nPts * realSize));
        spdlog::error("Storage requested exceeds device maximum variable size. Reason: {}. Try reducing storage to <{} time points, or reducing nPts. \n", nAux > nVar ? "aux vars" : "state vars", estimatedMaxStoreAlloc);
		throw std::invalid_argument("nPts*nStoreMax*nVar*realSize or nPts*nStoreMax*nAux*realSize is too big");
	}

	size_t currentTelements = currentStoreAlloc * nPts;

	//only resize device variables if size changed, or if not yet initialized
	if (!clInitialized || nStoreMax != currentStoreAlloc || telements != currentTelements)
	{

		nStoreMax = currentStoreAlloc;
		telements = currentTelements;
		xelements = nVar * currentTelements;
		auxelements = nAux * currentTelements;

		t.resize(telements);
		x.resize(xelements);
		dx.resize(xelements);
		aux.resize(auxelements);
		nStored.resize(nPts);

		//resize device variables
		try
		{
			//trajectory
			d_t = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * telements, NULL, &opencl.error);
			d_x = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * xelements, NULL, &opencl.error);
			d_dx = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * xelements, NULL, &opencl.error);
			d_aux = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, realSize * auxelements, NULL, &opencl.error);
			d_nStored = cl::Buffer(opencl.getContext(), CL_MEM_WRITE_ONLY, sizeof(int) * nPts, NULL, &opencl.error);
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODEtrajectory::resizeTrajectoryVariables():{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("resize d_t, d_x, d_dx, d_aux, d_nStored\n");
	}
}

//Simulation routine
void CLODEtrajectory::trajectory()
{
	if (clInitialized)
	{
		//resize output variables - will only occur if nPts or nSteps has changed [~4ms overhead on Tornado]
		resizeTrajectoryVariables();

		try
		{
			//kernel arguments
			int ix = 0;
			cl_trajectory.setArg(ix++, d_tspan);
			cl_trajectory.setArg(ix++, d_x0);
			cl_trajectory.setArg(ix++, d_pars);
			cl_trajectory.setArg(ix++, d_sp);
			cl_trajectory.setArg(ix++, d_xf);
			cl_trajectory.setArg(ix++, d_RNGstate);
			cl_trajectory.setArg(ix++, d_dt);
			cl_trajectory.setArg(ix++, d_t);
			cl_trajectory.setArg(ix++, d_x);
			cl_trajectory.setArg(ix++, d_dx);
			cl_trajectory.setArg(ix++, d_aux);
			cl_trajectory.setArg(ix++, d_nStored);

			//execute the kernel
			opencl.error = opencl.getQueue().enqueueNDRangeKernel(cl_trajectory, cl::NullRange, cl::NDRange(nPts));
			opencl.getQueue().finish();
		}
		catch (cl::Error &er)
		{
			spdlog::error("CLODEtrajectory::trajectory():{}({})\n", er.what(), CLErrorString(er.err()).c_str());
			throw er;
		}
		spdlog::debug("run trajectory\n");
	}
	else
	{
		printf("CLODE has not been initialized\n");
	}
}

std::vector<cl_double> CLODEtrajectory::getT()
{

	if (clSinglePrecision)
	{ //cast back to double
		std::vector<cl_float> tF(telements);
		opencl.error = copy(opencl.getQueue(), d_t, tF.begin(), tF.end());
		t.assign(tF.begin(), tF.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_t, t.begin(), t.end());
	}

	return t;
}

std::vector<cl_double> CLODEtrajectory::getX()
{

	if (clSinglePrecision)
	{ //cast back to double
		std::vector<cl_float> xF(xelements);
		opencl.error = copy(opencl.getQueue(), d_x, xF.begin(), xF.end());
		x.assign(xF.begin(), xF.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_x, x.begin(), x.end());
	}

	return x;
}

std::vector<cl_double> CLODEtrajectory::getDx()
{

	if (clSinglePrecision)
	{ //cast back to double
		std::vector<cl_float> dxF(xelements);
		opencl.error = copy(opencl.getQueue(), d_dx, dxF.begin(), dxF.end());
		dx.assign(dxF.begin(), dxF.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_dx, dx.begin(), dx.end());
	}

	return dx;
}

std::vector<cl_double> CLODEtrajectory::getAux()
{ //cast back to double

	if (clSinglePrecision)
	{
		std::vector<cl_float> auxF(auxelements);
		opencl.error = copy(opencl.getQueue(), d_aux, auxF.begin(), auxF.end());
		aux.assign(auxF.begin(), auxF.end());
	}
	else
	{
		opencl.error = copy(opencl.getQueue(), d_aux, aux.begin(), aux.end());
	}

	return aux;
}

std::vector<cl_int> CLODEtrajectory::getNstored()
{

	opencl.error = copy(opencl.getQueue(), d_nStored, nStored.begin(), nStored.end());
	return nStored;
}
