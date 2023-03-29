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
#include "observers.cl"

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
	prob.clRHSfilename="samples/lactotroph.cl";
	prob.nVar=8;
	prob.nPar=18;
	prob.nAux=1;
	prob.nWiener=1;
	prob.varNames.assign({"V", "n", "m", "b", "h", "h_T", "h_Na", "c"});
	prob.parNames.assign({"g_CaL", "g_CaT", "g_K", "g_SK", "g_Kir", "g_BK", "g_NaV", "g_A", "g_leak", "C_m", "E_leak","tau_m", "tau_ht", "tau_n", "tau_BK", "tau_h", "tau_hNa", "k_c"});
	prob.auxNames.assign({"i_noise"});
	
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

	//default pars
	std::vector<double> p({1.4, 0, 5, 0, 0, 0, 0, 0, 0.2, 10, -50, 1, 1, 39, 1, 1, 1, 0.03}); 
	
	// repeat the parameters nPts times: pack each paramater contiguously
	std::vector<double> pars(nPts, p[0]);
	for (int i=1; i<18; ++i)
		pars.insert(pars.end(), nPts, p[i]);

	//initial values: all zeros
	std::vector<double> x0(nPts*prob.nVar, 0.0);
	
	std::string observer = "thresh2";
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

	// Select device type and/or vendor using command line flags ("--device cpu/gpu/accel", "--vendor amd/intel/nvidia")
	OpenCLResource opencl( argc, argv);
	
	//prep timer and PRNG
	srand(static_cast <unsigned> (time(0))); 
    std::chrono::time_point<std::chrono::high_resolution_clock> start, end;
	std::chrono::duration<double, std::milli> elapsed_ms;
	
	
	// create the solver
	CLODEfeatures clo(prob, stepper, observer, CLSinglePrecision, opencl, CLODE_ROOT);
	
	clo.buildCL();
	
	//copy problem data to the device
	clo.initialize(tspan, x0, pars, sp, op); 
	
	// clo.seedRNG(mySeed);
	
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
