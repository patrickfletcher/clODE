
#include "clODE_random.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"
#include "steppers.cl"

__kernel void initializeObserver(
	__constant realtype *tspan,			//time vector [t0,tf] - adds (tf-t0) to these at the end
	__global realtype *x0,				//initial state 				[nPts*nVar]
	__constant realtype *pars,			//parameter values				[nPts*nPar]
	__constant struct SolverParams *sp, //dtmin/max, tols, etc
	__global ulong *RNGstate,			//enables host seeding/continued streams	    [nPts*nRNGstate]
    __global realtype *d_dt,            //array of dt values, one per solver
	__global ObserverData *OData,		//for continue
	__constant struct ObserverParams *opars)
{
	int i = get_global_id(0);
	int nPts = get_global_size(0);

	realtype ti, dt;
    realtype p[N_PAR], xi[N_VAR], dxi[N_VAR], auxi[N_AUX], wi[N_WIENER];
	rngData rd;

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

    for (int j = 0; j < N_WIENER; ++j)
#ifdef STOCHASTIC_STEPPER
        wi[j] = randn(&rd) / sqrt(dt);
#else
        wi[j] = RCONST(0.0);
#endif
	getRHS(ti, xi, p, dxi, auxi, wi); //slope at initial point, needed for FSAL

	ObserverData odata = OData[i]; //private copy of observer data

	initializeObserverData(&ti, xi, dxi, auxi, &odata, opars);

#ifdef TWO_PASS_EVENT_DETECTOR

    int step = 0;
    int stepflag = 0;
	while (ti < tspan[1] && step < sp->max_steps)
	{
		++step;
        stepflag = stepper(&ti, xi, dxi, p, sp, &dt, tspan, auxi, wi, &rd);
        // if (stepflag!=0)
        //     break;

		warmupObserverData(&ti, xi, dxi, auxi, &odata, opars);
	}

#endif //TWO_PASS_EVENT_DETECTOR

	initializeEventDetector(&ti, xi, dxi, auxi, &odata, opars);

	//update the global ObserverData array
	OData[i] = odata;

	//dt only evolves for TWO_PASS_EVENT_DETECTOR, in which case we want to restart. Don't save dt.
}
