#include "clODE_random.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "observers.cl"
#include "steppers.cl"
#include "realtype.cl"

__kernel void initializeObserver(
	__constant realtype *tspan,			//time vector [t0,tf] - adds (tf-t0) to these at the end
	__global realtype *x0,				//initial state 				[nPts*nVar]
	__constant realtype *pars,			//parameter values				[nPts*nPar]
	__constant struct SolverParams *sp, //dtmin/max, tols, etc
	__global ulong *RNGstate,			//enables host seeding/continued streams	    [nPts*nRNGstate]
	__global ObserverData *OData,		//for continue
	__constant struct ObserverParams *opars)
{

	int i = get_global_id(0);
	int nPts = get_global_size(0);

	realtype ti, dt;
	realtype xi[N_VAR], dxi[N_VAR], auxi[N_AUX];
	realtype p[N_PAR], wi[N_WIENER];
	rngData rd;
	int step, stepflag;
	__constant realtype *tspanPtr = tspan;

	//get private copy of ODE parameters, initial data, and compute slope at initial state
	ti = tspan[0];
	dt = sp->dt;

	for (int j = 0; j < N_PAR; ++j)
		p[j] = pars[j * nPts + i];

	for (int j = 0; j < N_VAR; ++j)
		xi[j] = x0[j * nPts + i];

	for (int j = 0; j < N_RNGSTATE; ++j)
		rd.state[j] = RNGstate[j * nPts + i];

	rd.randnUselast = 0;

#ifdef ADAPTIVE_STEPSIZE
	realtype lastErr=RCONST(1.0);
	realtype lastDtRatio=RCONST(1.0);
	for (int j = 0; j < N_WIENER; ++j)
		wi[j] = RCONST(0.0);
#else
	for (int j = 0; j < N_WIENER; ++j)
		wi[j] = randn(&rd) / sqrt(dt);
#endif

	getRHS(ti, xi, p, dxi, auxi, wi); //slope at initial point, needed for FSAL steppers (bs23, dorpri5)

	ObserverData odata = OData[i]; //private copy of observer data


	initializeObserverData(&ti, xi, dxi, auxi, &odata, opars);

#ifdef TWO_PASS_EVENT_DETECTOR

	step = 0;
	stepflag = 0;
	while (ti < tspan[1] && step < sp->max_steps)
	{
		++step;
# ifdef ADAPTIVE_STEPSIZE
		//leave the wi=0 for adaptive steppers
		stepflag = stepper(&ti, xi, dxi, p, sp, &dt, tspanPtr, auxi, wi, &lastErr, &lastDtRatio);
		if (stepflag!=0)
			break;
# else
		//update Wiener variables - fixed size steppers can scale by dt here
		for (int j = 0; j < N_WIENER; ++j)
			wi[j] = randn(&rd) / sqrt(dt); //NOTE: divide by sqrt(dt) because Euler will multiply this by dt in the stepper.

		stepper(&ti, xi, dxi, p, dt, auxi, wi);
		ti = tspan[0] + step * dt; //purify ti - Gets nSteps correct, but incompatible with shrinking final step without conditional to check if doing the last step
# endif

		//FSAL: dxi is at new ti, Not FSAL: dxi is at old ti
		warmupObserverData(&ti, xi, dxi, auxi, &odata, opars);
	}

#endif //TWO_PASS_EVENT_DETECTOR

	initializeEventDetector(&ti, xi, dxi, auxi, &odata, opars);

	//update the global ObserverData array
	OData[i] = odata;
}
