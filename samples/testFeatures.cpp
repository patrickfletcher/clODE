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

#define CLODE_ROOT "clode/cpp/"

//Generate random points within given bounds
template<typename T> std::vector<T> generateRandomPoints(std::vector<T> lb, std::vector<T> ub, int nPts);

//currently the only command line arguments are to select device/vendor type ("--device cpu/gpu/accel", "--vendor amd/intel/nvidia")
int main(int argc, char **argv)
{
	try 
	{
 	cl_int nPts=32;
	bool CLSinglePrecision=true;
	
	ProblemInfo prob;
    prob.clRHSfilename="samples/lactotroph.cl";;
	prob.nVar=4;
	prob.nPar=3;
	prob.nAux=1;
	prob.nWiener=0;
	prob.varNames.assign({"v","n","f","c"});
	prob.parNames.assign({"gcal","gsk","gbk"});
	prob.auxNames.assign({"ical"});
	
	std::string stepper="rk4";
	
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

	int mySeed=1;

	//default pars
	std::vector<double> p({1.5,3.0,1.0});
	
	// repeat the parameters nPts times: pack each paramater contiguously
	std::vector<double> pars(nPts, p[0]);
	pars.insert(pars.end(), nPts, p[1]);
	pars.insert(pars.end(), nPts, p[2]);

	// //Parameter sets will be sampled uniformly from [lb,ub] for each parameter
	// std::vector<double> lb({1.0,5.0,1.0});
	// std::vector<double> ub({1.0,5.0,1.0});
	// std::vector<double> pars=generateRandomPoints(lb, ub, nPts);
	
	//initial values: all zeros
	std::vector<double> x0(nPts*prob.nVar, 0.0);
	
	std::string observer = "thresh2";
	ObserverParams<double> op;
	op.eVarIx=0;
	op.fVarIx=0;
	op.maxEventCount=100;
	op.maxEventTimestamps=0;
	op.minXamp=0.0;
	op.minIMI=0.0;
	op.nHoodRadius=0.01;
	op.xUpThresh=0.3;
	op.xDownThresh=0.2;
	op.dxUpThresh=0.0;
	op.dxDownThresh=0.0;
	op.eps_dx=1e-7;

	
	// Select device type and/or vendor using command line flags ("--device cpu/gpu/accel", "--vendor amd/intel/nvidia")
	OpenCLResource opencl( argc, argv);
	
	//prep timer and PRNG
	srand(static_cast <unsigned> (time(0))); 
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double, std::milli> elapsed_ms;

	// create the simulator
	CLODEfeatures clo(prob, stepper, observer, op, CLSinglePrecision, opencl, CLODE_ROOT);

	//~ std::cout<<"here"<<std::endl;
	//copy problem data to the device

    clo.buildCL();
	clo.setSolverParams(sp);
	clo.setProblemData(x0, pars);
	clo.setTspan(tspan);
	clo.seedRNG(mySeed);
	
	//run the simulation 
	clo.transient();
    clo.shiftX0();


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
	
	std::vector<std::string> fnames=clo.getFeatureNames();

	std::cout<< "\ntf="<< tspan[0] <<std::endl;
	std::cout<<"Features:"<< "\n";
	for (int i=0; i<clo.getNFeatures(); ++i)
		std::cout << " " << fnames[i] << " = " << F[i*nPts] << "\n";
	
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
