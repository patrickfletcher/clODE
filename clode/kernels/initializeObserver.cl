
#include "clODE_random.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"
#include "steppers.cl"

__kernel void initializeObserver(
    __constant realtype *tspan,         //time interval
    __global realtype *x0,              //initial state 	   [nPts*nVar]
    __constant realtype *pars,          //parameter values	   [nPts*nPar]
    __constant struct IntegrationSettings *settings, //dtmin/max, tols
    __global ulong *RNGstate,           //final RNG	state	   [nPts*nRNGstate]
	__global realtype *RNGspareNormal,  //cached Box-Muller spare normal [nPts]
	__global uint *RNGspareNormalValid, //cached Box-Muller availability [nPts]
	__global realtype *preparedWiener,  //prepared next-step Wiener sample [nPts*nWiener]
	__global uint *preparedWienerValid, //prepared next-step Wiener availability [nPts]
    __global realtype *d_dt,            //final dt values      [nPts]
	__global ObserverState *observer_states, //persistent observer state
	__constant struct ObserverRuntimeSettings *opars) //observer runtime settings
{
	int i = get_global_id(0);
	int nPts = get_global_size(0);

	realtype ti, dt, solveElapsed, solveElapsedCorrection;
    realtype p[N_PAR], xi[N_VAR], dxi[N_VAR];
    realtype auxi[N_AUX>0?N_AUX:1];
    realtype wi[N_WIENER>0?N_WIENER:1];
	struct rngData rd;

	//get private copy of ODE parameters, initial data, and compute slope at initial state
	ti = tspan[0];
    dt = d_dt[i];
	solveElapsed = ZERO;
	solveElapsedCorrection = ZERO;
	realtype solveDuration = tspan[1] - tspan[0];

	for (int j = 0; j < N_PAR; ++j)
		p[j] = pars[j * nPts + i];

	for (int j = 0; j < N_VAR; ++j)
		xi[j] = x0[j * nPts + i];

	for (int j = 0; j < N_RNGSTATE; ++j)
		rd.state[j] = RNGstate[j * nPts + i];

	rd.randnUselast = RNGspareNormalValid[i] != 0;
	rd.randnLast = RNGspareNormal[i];

	ObserverState observer_state = observer_states[i];

    // generate random numbers if needed
    for (int j = 0; j < N_WIENER; ++j)
#ifdef STOCHASTIC_STEPPER
		wi[j] = preparedWienerValid[i] != 0 ? preparedWiener[j * nPts + i] : randn(&rd) / sqrt(dt);
#else
        wi[j] = ZERO;
#endif

    //get the slope and aux at initial point
    getRHS(ti, xi, p, dxi, auxi, wi); 

	initializeObserverState(&ti, xi, dxi, auxi, &observer_state, opars);

#ifdef TWO_PASS_EVENT_DETECTOR

	//time-stepping loop
	ulong step = 0;
	realtype acceptedStepDt = ZERO;
	while (
		compensatedTimeValue(solveElapsed, solveElapsedCorrection) < solveDuration
		&& step < settings->max_steps
	)
	{
		++step;
		int stepflag = stepper(
			&ti,
			&solveElapsed,
			&solveElapsedCorrection,
			xi,
			dxi,
			p,
			settings,
			&dt,
			&acceptedStepDt,
			tspan,
			auxi,
			wi,
			&rd
		);
		if (stepflag != 0)
            break;

		warmupObserverState(&ti, xi, dxi, auxi, &observer_state, opars);
	}
	//rewind the time and state so initializeEventDetector gets the right values
	ti = tspan[0];
	solveElapsed = ZERO;
	solveElapsedCorrection = ZERO;
	for (int j = 0; j < N_VAR; ++j)
		xi[j] = x0[j * nPts + i];
	getRHS(ti, xi, p, dxi, auxi, wi);

#endif //TWO_PASS_EVENT_DETECTOR

	initializeEventDetector(&ti, xi, dxi, auxi, &observer_state, opars);

	// Update the persistent observer state array.
	observer_states[i] = observer_state;

}
