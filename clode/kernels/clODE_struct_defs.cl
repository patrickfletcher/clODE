#ifndef CLODE_STRUCT_DEFS_H_
#define CLODE_STRUCT_DEFS_H_

//TODO: bounds param (like XPP) to catch numerical instability/blow up in finite time?
//TODO: provide different structs for base vs trajectory solvers (expose relevant members)
#include "realtype.cl"

struct IntegrationSettings
{
	realtype dt;
	realtype dtmax;
	realtype abstol;
	realtype reltol;
	ulong max_steps;
};

struct TrajectoryOutputSettings
{
	unsigned int max_store;
	unsigned int nout;
};

struct ObserverRuntimeSettings
{
	unsigned int eVarIx; //variable for event detection
	unsigned int fVarIx; //variable for features

	unsigned int maxEventCount; //time loop limiter
	unsigned int eventDirection; //0=rising, 1=falling, 2=either
	// maxEventTimestamps - not here: used as #define N_STORE_EVENTS for fixed-size event timestamp buffer
	realtype minXamp;  //consider oscillations lower than this to be steady state (return mean X)
	realtype minIMI;

	//neighborhood return map
	realtype nHoodRadius;

	//section. Two interpretations: absolute, relative
	realtype xUpThresh;
	realtype xDownThresh;
	realtype dxUpThresh;
	realtype dxDownThresh;

	//local extremum - tolerance for zero crossing of dx - for single precision: if RHS involves sum of terms of O(1), dx=zero is noise at O(1e-7)
	realtype eps_dx;
};

#endif //CLODE_STRUCT_DEFS_H_
