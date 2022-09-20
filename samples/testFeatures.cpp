/*
 * main.cpp: test example to run clODE and clODEtrajectory.cpp
 * 
 * Copyright 2017 Patrick Fletcher <patrick.fletcher@nih.gov>
 * 
 */
 
 //TODO what is the right way to do unit testing?
 
#include <algorithm>
#include <cfloat>
#include <chrono>
#include <iostream>
#include <fstream>

#include "OpenCLResource.hpp"
#include "CLODE.hpp"
#include "CLODEfeatures.hpp"

#define CLODE_ROOT "src/"

//Generate random points within given bounds
template<typename T> std::vector<T> generateRandomPoints(std::vector<T> lb, std::vector<T> ub, int nPts);

//currently the only command line arguments are to select device/vendor type ("--device cpu/gpu/accel", "--vendor amd/intel/nvidia")
int main(int argc, char **argv)
{
 	#if defined(WIN32)||defined(_WIN64)
 		_putenv_s("CUDA_CACHE_DISABLE", "1");
 	#else
 		setenv("CUDA_CACHE_DISABLE", "1", 1);
 	#endif
	try 
	{
		
	cl_int nPts=4096;
	bool CLSinglePrecision=true;
	
	ProblemInfo prob;

	prob.clRHSfilename="samples/lactotroph.cl";
		
	prob.nVar=4;
	prob.nPar=3;
	prob.nAux=1;
	prob.nWiener=1;

    prob.varNames.push_back("a");
    prob.varNames.push_back("b");
    prob.varNames.push_back("c");
    prob.varNames.push_back("d");

    prob.parNames.push_back("aa");
    prob.parNames.push_back("bb");
    prob.parNames.push_back("cc");

    prob.auxNames.push_back("dd");


	std::string stepper="euler";
	
	//parameters for solver and objective function
	std::vector<double> tspan({0.0,1000.0});
	
	int nReps=1;
	
	SolverParams<double> sp;
	sp.dt=0.1;
	sp.dtmax=1.00;
	sp.abstol=1e-6;
	sp.reltol=1e-3;
	sp.max_steps=10000000;
	sp.max_store=10000000;
	sp.nout=50;
	
	std::string observer="basic";
	
	ObserverParams<double> op;
	op.eVarIx=0;
	op.fVarIx=0;
	op.maxEventCount=100;
	op.minXamp=1;
	op.nHoodRadius=0.01;
	op.xUpThresh=0.3;
	op.xDownThresh=0.2;
	op.dxUpThresh=0;
	op.dxDownThresh=0;
	op.eps_dx=1e-7;
	
	
	int mySeed=1;
	
	//default pars
	std::vector<double> p({1.0,5.0,1.0}); 
	
	//Parameter sets will be sampled uniformly from [lb,ub] for each parameter
	std::vector<double> lb({1.0,5.0,1.0});
	std::vector<double> ub({1.0,5.0,1.0});
	std::vector<double> pars=generateRandomPoints(lb, ub, nPts);
	
	//initial values: all zeros
	std::vector<double> x0(nPts*prob.nVar, 0.0);

	
	// Select device type and/or vendor using command line flags ("--device cpu/gpu/accel", "--vendor amd/intel/nvidia")
	OpenCLResource opencl( argc, argv);
	
	//prep timer and PRNG
	srand(static_cast <unsigned> (time(0))); 
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double, std::milli> elapsed_ms;


	// create the simulator
	CLODEfeatures clo(prob, stepper, observer, CLSinglePrecision, opencl, CLODE_ROOT);

	//~ std::cout<<"here"<<std::endl;
	//copy problem data to the device

	clo.initialize(tspan, x0, pars, sp, op);

	clo.seedRNG(mySeed);
	
	//run the simulation 
	clo.transient();

		
	start = std::chrono::high_resolution_clock::now();
		
	std::cout<<std::endl;
	for(int i=0;i<nReps;++i){
		clo.features();
	}
	
	end = std::chrono::high_resolution_clock::now();
	elapsed_ms += end-start;
	
	//retrieve result from device
	std::vector<double> F=clo.getF();
	tspan=clo.getTspan();
	std::cout<< "\ntf="<< tspan[0] <<", F:"<< "\n";
	for (int i=0; i<clo.getNFeatures(); ++i)
		std::cout<< " " << F[i*nPts];
	
	std::cout<<std::endl;
	
    std::cout<< "Compute time: " << elapsed_ms.count() << "ms\n";
	std::cout<<std::endl;
	
	
	} catch (std::exception &er) {
        std::cout<< "ERROR: " << er.what() << std::endl;
        std::cout<<"exiting...\n";
		return -1;
	}
    
	return 0;
}


//Generate random points within given bounds. Pack coordinates contiguously: all x1, then x2, etc.
template<typename T> std::vector<T> generateRandomPoints(std::vector<T> lb, std::vector<T> ub, int nPts)
{
	int dim=lb.size();
	std::vector<T> x(nPts*dim);
	
	T r;
	for (int i=0; i<dim; ++i)
	{	
		for (int j=0; j<nPts; ++j)
		{
			r = static_cast <T> (rand()) / static_cast <T> (RAND_MAX); //in [0,1]
			x[i*nPts+j]=lb[i] +r*(ub[i]-lb[i]); 
		}
	}
	return x;
}

//explicit instantiation of template function
template std::vector<float> generateRandomPoints(std::vector<float> lb, std::vector<float> ub, int nPts);
template std::vector<double> generateRandomPoints(std::vector<double> lb, std::vector<double> ub, int nPts);
