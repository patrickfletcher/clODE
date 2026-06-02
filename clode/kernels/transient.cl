
#include "clODE_random.cl"
#include "clODE_solver_status.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "realtype.cl"
#include "steppers.cl"

// the most basic trajectory solver that stores nothing but the final variable values (and RNG state)
__kernel void transient(
    __constant realtype *tspan,         //time interval
    __global realtype *t0,              //per-item absolute chunk start [nPts]
    __global realtype *x0,              //initial state 	   [nPts*nVar]
    __constant realtype *pars,          //parameter values	   [nPts*nPar]
    __constant struct IntegrationSettings *settings, //dtmin/max, tols
    __global realtype *xf,              //final state 		   [nPts*nVar]
    __global ulong *RNGstate,           //final RNG	state	   [nPts*nRNGstate]
    __global realtype *RNGspareNormal,  //cached Box-Muller spare normal [nPts]
    __global uint *RNGspareNormalValid, //cached Box-Muller availability [nPts]
    __global realtype *preparedWiener,  //prepared next-step Wiener sample [nPts*nWiener]
    __global uint *preparedWienerValid, //prepared next-step Wiener availability [nPts]
    __global realtype *d_dt,            //final dt values      [nPts]
    __global int *status,               //solve status values  [nPts]
    __global ulong *stepCount,          //accepted step counts [nPts]
    __global realtype *acceptedDt,      //last accepted dt     [nPts]
    __global realtype *tf)              //final time values    [nPts]
{
    int i = get_global_id(0);
    int nPts = get_global_size(0);

    realtype ti, dt, solveElapsed, solveElapsedCorrection;
    realtype p[N_PAR], xi[N_VAR], dxi[N_VAR];
    realtype auxi[N_AUX>0?N_AUX:1];
    realtype wi[N_WIENER>0?N_WIENER:1];
    struct rngData rd;

    //get private copy of ODE parameters, initial data, and compute slope at initial state
    realtype tOrigin = t0[i];
    ti = tOrigin;
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

    // generate random numbers if needed
    for (int j = 0; j < N_WIENER; ++j)
#ifdef STOCHASTIC_STEPPER
        wi[j] = preparedWienerValid[i] != 0 ? preparedWiener[j * nPts + i] : randn(&rd) / sqrt(dt);
#else
        wi[j] = ZERO;
#endif

    //get the slope and aux at initial point
    getRHS(ti, xi, p, dxi, auxi, wi); 

	//time-stepping loop
    ulong step = 0;
    ulong acceptedSteps = 0;
    realtype acceptedStepDt = ZERO;
    realtype lastAcceptedStepDt = ZERO;
    int solveStatus = SOLVE_STATUS_COMPLETED;
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
            tOrigin,
            solveDuration,
            auxi,
            wi,
            &rd
        );
        if (stepflag != 0)
        {
            solveStatus = (
                stepflag == SOLVE_STATUS_NO_PROGRESS
                ? SOLVE_STATUS_NO_PROGRESS
                : SOLVE_STATUS_STEPPER_FAILED
            );
            break;
        }
        acceptedSteps = step;
        lastAcceptedStepDt = acceptedStepDt;
    }

    solveStatus = finalizeSolveStatus(
        solveStatus,
        compensatedTimeValue(solveElapsed, solveElapsedCorrection) >= solveDuration,
        step >= settings->max_steps,
        false,
        false
    );

    //write the final solution values to global memory.
    for (int j = 0; j < N_VAR; ++j)
        xf[j * nPts + i] = xi[j];

    // To get same RNG on repeat (non-continued) run, need to set the seed to same value
    for (int j = 0; j < N_RNGSTATE; ++j)
        RNGstate[j * nPts + i] = rd.state[j];

    RNGspareNormal[i] = rd.randnLast;
    RNGspareNormalValid[i] = rd.randnUselast ? 1 : 0;

#ifdef STOCHASTIC_STEPPER
    for (int j = 0; j < N_WIENER; ++j)
        preparedWiener[j * nPts + i] = wi[j];
    preparedWienerValid[i] = 1;
#else
    preparedWienerValid[i] = 0;
#endif

    // update dt to its final value (for adaptive stepper continue)
    d_dt[i] = dt;

    status[i] = solveStatus;
    stepCount[i] = acceptedSteps;
    acceptedDt[i] = lastAcceptedStepDt;

    // store the actual final time value
    tf[i] = ti;
}
