//TODO: investigate sources of roundoff error in long "continues" - running means?

#include "clODE_utilities.h"
#include "clODE_random.h"
#include "clODE_struct_defs.h"
#include "observers.h"
#include "realtype.h"
#include "steppers.h"

__kernel void features(
__constant   realtype * tspan,			//time vector [t0,tf] - adds (tf-t0) to these at the end	
__global     realtype * x0,				//initial state 				[nPts*nVar]
__constant   realtype * pars,			//parameter values				[nPts*nPar]
__constant   struct SolverParams * sp,	//dtmin/max, tols, etc
__global     ulong * RNGstate,          //enables host seeding/continued streams	    [nPts*nRNGstate]
__global   	 ObserverData * OData,	    //for continued feature detection
__constant   struct ObserverParams * opars,
int doInitialization,               //reset observerdata?
__global     realtype * xf,				//final state 				[nPts*nVar]
__global     realtype * F,          //feature results
__global     realtype * t,			// trajectory storage: t, x, dx, aux
__global     realtype * x,			
__global     realtype * dx,			
__global     realtype * aux,		
			 int nStoreMax,
__global	 int * nStored
) {

int i=get_global_id(0); 
int nPts=get_global_size(0); 

realtype ti, dt; 
realtype xi[N_VAR],dxi[N_VAR],auxi[N_AUX];
realtype p[N_PAR], wi[N_WIENER];
rngData rd;
int step, stepflag;

//get private copy of ODE parameters, initial data, and compute slope at initial state
ti = tspan[0];
dt=sp->dt;

for (int j=0; j<N_PAR;++j) 
	p[j]=pars[j*nPts + i];

for (int j=0; j<N_VAR;++j) 
	xi[j] = x0[j*nPts + i];

for (int j=0; j<N_RNGSTATE; ++j)
	rd.state[j]=RNGstate[j*nPts+i];

rd.randnUselast=0;

#ifdef ADAPTIVE_STEPSIZE 
for (int j=0; j<N_WIENER; ++j)
	wi[j]=RCONST(0.0); 
#else
for (int j=0; j<N_WIENER; ++j)
	wi[j]=randn(&rd)/sqrt(dt);   
#endif    
	
getRHS(ti, xi, p, dxi, auxi, wi); //slope at initial point, needed for FSAL steppers (bs23, dorpri5)
 
//store the initial point

int storeix=0;
t[storeix+i]=tspan[0];
for (int j=0; j<N_VAR;++j) 
    x[storeix*nPts*N_VAR+j*nPts+i]=xi[j];

for (int j=0; j<N_VAR;++j) 
    dx[storeix*nPts*N_VAR+j*nPts+i]=dxi[j];

for (int j=0; j<N_AUX;++j) 
    aux[storeix*nPts*N_AUX+j*nPts+i]=auxi[j];
    

ObserverData odata = OData[i]; //private copy of observer data


if (doInitialization==1) {
	
	initializeObserverData(&ti, xi, dxi, auxi, &odata, opars);
	
#ifdef TWO_PASS_EVENT_DETECTOR
	
	step=0;
	stepflag=0;
	while (ti < tspan[1] && step<sp->max_steps && stepflag==0) {
		
		++step;
	#ifdef ADAPTIVE_STEPSIZE
		//leave the wi=0 for adaptive steppers
		stepflag=stepper(&ti, xi, dxi, p, sp, &dt, tspan, auxi, wi); 
	#else
		//update Wiener variables - fixed size steppers can scale by dt here
		for (int j=0; j<N_WIENER; ++j)
			wi[j]=randn(&rd)/sqrt(dt);   //NOTE: divide by sqrt(dt) because Euler will multiply this by dt in the stepper.
			
		stepper(&ti, xi, dxi, p, dt, auxi, wi);
		ti=tspan[0] + step*dt;   //purify ti - Gets nSteps correct, but incompatible with shrinking final step without conditional to check if doing the last step
	#endif
		
		//FSAL: dxi is at new ti, Not FSAL: dxi is at old ti
		warmupObserverData(&ti, xi, dxi, auxi, &odata, opars);
	}
	
	//rewind state to t0, x0, RNG0
	//or just reset ti->t0?
	ti = tspan[0];

	// dt=sp->dt;

	// for (int j=0; j<N_VAR;++j) 
	// 	xi[j] = x0[j*nPts + i];

	// for (int j=0; j<N_RNGSTATE; ++j)
	// 	rd.state[j]=RNGstate[j*nPts+i];

	// rd.randnUselast=0;
		
	// #ifdef ADAPTIVE_STEPSIZE 
	// for (int j=0; j<N_WIENER; ++j)
	// 	wi[j]=RCONST(0.0);  
	// #else
	// for (int j=0; j<N_WIENER; ++j)
	// 	wi[j]=randn(&rd)/sqrt(dt);   
	// #endif   
		
	// getRHS(ti, xi, p, dxi, auxi, wi);
	
#endif //TWO_PASS_EVENT_DETECTOR

	initializeEventDetector(&ti, xi, dxi, auxi, &odata, opars);
}
//else { //no initialization: prepare observer data for continuation if needed
	//prepareObserverData(&ti, xi, dxi, auxi, &odata, opars);
//}



//time-stepping loop, main time interval
step=0;
stepflag=0;
bool eventOccurred;
bool terminalEvent;
while (ti < tspan[1] && step<sp->max_steps && stepflag==0) {
	
	++step;
	
#ifdef ADAPTIVE_STEPSIZE
	//leave the wi=0 for adaptive steppers
	stepflag=stepper(&ti, xi, dxi, p, sp, &dt, tspan, auxi, wi); 
#else
	//update Wiener variables - fixed size steppers scale by dt here
	for (int j=0; j<N_WIENER; ++j)
		wi[j]=randn(&rd)/sqrt(dt);   //NOTE: divide by sqrt(dt) because Euler will multiply this by dt in the stepper.
		
	stepper(&ti, xi, dxi, p, dt, auxi, wi);
	ti=tspan[0] + step*dt;   //purify ti - Gets nSteps correct
#endif

	++odata.stepcount;

	//FSAL: dxi is at new ti, Not FSAL: dxi is at old ti
	eventOccurred=eventFunction(&ti, xi, dxi, auxi, &odata, opars);
	if (eventOccurred) {
		terminalEvent=computeEventFeatures(&ti, xi, dxi, auxi, &odata, opars);
		if (terminalEvent){break;};
	}
	
	updateObserverData(&ti, xi, dxi, auxi, &odata, opars, eventOccurred);  //TODO: if not FSAL, dxi buffer is delayed by one. (dxi is slope at LAST timestep)

    
	//store every sp.nout'th step after the initial point
	if (step%sp->nout==0) 
	{
		++storeix;
		
		t[storeix*nPts+i]=ti; //adaptive steppers give different timepoints for each trajectory
		
		for (int j=0; j<N_VAR;++j)
			x[storeix*nPts*N_VAR+j*nPts+i]=xi[j];
        
		for (int j=0; j<N_VAR;++j)
#ifdef FSAL_STEP_PROPERTY
			dx[storeix*nPts*N_VAR+j*nPts+i]=dxi[j];
#else
			dx[(storeix-1)*nPts*N_VAR+j*nPts+i]=dxi[j]; //without FSAL, dxi returned by stepper is evaluated at beginning of timestep
#endif
		
		for (int j=0; j<N_AUX;++j)
			aux[storeix*nPts*N_AUX+j*nPts+i]=auxi[j];
		
	}
}

//readout features of interest and write to global F:
finalizeFeatures(&ti, xi, dxi, auxi, &odata, opars, F, i, nPts);

//finalize observerdata for possible continuation
finalizeObserverData(&ti, xi, dxi, auxi, &odata, opars, tspan);


#if !defined(FSAL_STEP_PROPERTY)
//get the slope at the last point
for (int j=0; j<N_WIENER; ++j)
	wi[j]=randn(&rd)/sqrt(dt);

getRHS(ti, xi, p, dxi, auxi, wi);

for (int j=0; j<N_VAR;++j)
	dx[storeix*nPts*N_VAR+j*nPts+i]=dxi[j];
#endif

nStored[i]=storeix; //storeix ranged from 0 to nStored-1

//set global variables to be ready to continue
for (int j=0; j<N_VAR;++j)
	xf[j*nPts + i]=xi[j];

for (int j=0; j<N_RNGSTATE; ++j)
	RNGstate[j*nPts+i]=rd.state[j];

//TODO: odataf - store separate observerdata final state to match xf shift model???
OData[i] = odata;
}