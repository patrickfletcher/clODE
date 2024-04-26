// trajectory solver that enables observers

#include "clODE_random.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"
#include "steppers.cl"

__kernel void features(
	__constant realtype *tspan,         //time vector [t0,tf] - adds (tf-t0) to these at the end
	__global realtype *x0,              //initial state 				[nPts*nVar]
	__constant realtype *pars,          //parameter values				[nPts*nPar]
	__constant struct SolverParams *sp, //dtmin/max, tols, etc
	__global realtype *xf,              //final state 				[nPts*nVar]
	__global ulong *RNGstate,           //state for RNG					[nPts*nRNGstate]
    __global realtype *d_dt,            //array of dt values, one per solver
	__global ObserverData *OData,		//Observer data
	__constant struct ObserverParams *opars,
	__global realtype *F)
{
	int i = get_global_id(0);
	int nPts = get_global_size(0);

	realtype ti, dt;
    realtype p[N_PAR], xi[N_VAR], dxi[N_VAR];
    realtype auxi[N_AUX>0?N_AUX:1];
    realtype wi[N_WIENER>0?N_WIENER:1];
	struct rngData rd;

	//get private copy of ODE parameters, initial data, and compute slope at initial state
	ti = tspan[0];
    dt = d_dt[i];

	for (int j = 0; j < N_PAR; ++j)
		p[j] = pars[j * nPts + i];

	for (int j = 0; j < N_VAR; ++j)
		xi[j] = x0[j * nPts + i];

	for (int j = 0; j < N_RNGSTATE; ++j)
		rd.state[j] = RNGstate[j * nPts + i];

	ObserverData odata = OData[i]; 

    // generate random numbers if needed
    rd.randnUselast = 0;
    for (int j = 0; j < N_WIENER; ++j)
#ifdef STOCHASTIC_STEPPER
        wi[j] = randn(&rd) / sqrt(dt);
#else
        wi[j] = ZERO;
#endif

    //get the slope and aux at initial point
    getRHS(ti, xi, p, dxi, auxi, wi); 

	//time-stepping loop
    unsigned int step = 0;
    int stepflag = 0;
	bool eventOccurred;
	bool terminalEvent;
	while (ti <= tspan[1] && step < sp->max_steps)
	{
		++step;
        stepflag = stepper(&ti, xi, dxi, p, sp, &dt, tspan, auxi, wi, &rd);
        // if (stepflag!=0)
            // break;

		updateObserverData(&ti, xi, dxi, auxi, &odata, opars); 

		eventOccurred = eventFunction(&ti, xi, dxi, auxi, &odata, opars);
		if (eventOccurred)
		{
			// record current state into dedicated event buffers (same pattern as trajectory)

			terminalEvent = computeEventFeatures(&ti, xi, dxi, auxi, &odata, opars);
			if (terminalEvent)
				break;
		}
	}

	//readout features of interest and write to global F:
	finalizeFeatures(&ti, xi, dxi, auxi, &odata, opars, F, i, nPts);

	//finalize observerdata for possible continuation
	finalizeObserverData(&ti, xi, dxi, auxi, &odata, opars, tspan);

	//store the observerData in global memory
	OData[i] = odata;

    //write the final solution values to global memory.
	for (int j = 0; j < N_VAR; ++j)
		xf[j * nPts + i] = xi[j];

    // To get same RNG on repeat (non-continued) run, need to set the seed to same value
	for (int j = 0; j < N_RNGSTATE; ++j)
		RNGstate[j * nPts + i] = rd.state[j];

    // update dt to its final value (for adaptive stepper continue)
    d_dt[i] = dt;
}
