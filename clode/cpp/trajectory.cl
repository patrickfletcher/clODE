//trajectory: stores in global variables directly

//TODO: optionally store only a subset of variables (allow bigger nPts): sp.varIx vector of indices? already in observerpars...???
//TODO: alternate storage at specified time points only - host sets t vector, interp and store x/dx/aux whenever ti passes t[nextstoreix]
//TODO: support dense output (refine option) for solvers like Dopri45 (or all if we force FSAL; use interpolant of same order as solver)
//TODO: is there any way to avoid writing to global at each store step? shared mem?

#include "clODE_random.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "realtype.cl"
#include "steppers.cl"

__kernel void trajectory(
    __constant realtype *tspan,         //time vector [t0,tf] - adds (tf-t0) to these at the end
    __global realtype *x0,              //initial state 				[nPts*nVar]
    __constant realtype *pars,          //parameter values				[nPts*nPar]
    __constant struct SolverParams *sp, //dtmin/max, tols, etc
    __global realtype *xf,              //final state 				[nPts*nVar]
    __global ulong *RNGstate,           //state for RNG				[nPts*nRNGstate]
    __global realtype *d_dt,            //array of dt values, one per solver
    __global realtype *t,               //
    __global realtype *x,               //
    __global realtype *dx,              //
    __global realtype *aux,             //
    __global int *nStored)
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

    //store the initial point
    int storeix = 0;
    t[storeix * nPts + i] = ti;
    for (int j = 0; j < N_VAR; ++j)
        x[storeix * nPts * N_VAR + j * nPts + i] = xi[j];
    for (int j = 0; j < N_VAR; ++j)
        dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];
    for (int j = 0; j < N_AUX; ++j)
        aux[storeix * nPts * N_AUX + j * nPts + i] = auxi[j];
    
	//time-stepping loop
    unsigned int step = 0;
    int stepflag = 0;
    while (ti <= tspan[1] && step < sp->max_steps && storeix < sp->max_store)
    {
		++step;
        stepflag = stepper(&ti, xi, dxi, p, sp, &dt, tspan, auxi, wi, &rd);
        // if (stepflag!=0)
        //     break;

        //store every sp.nout'th step after the initial point
        if (step % sp->nout == 0)
        {
            ++storeix;
            t[storeix * nPts + i] = ti;
            for (int j = 0; j < N_VAR; ++j)
                x[storeix * nPts * N_VAR + j * nPts + i] = xi[j];
            for (int j = 0; j < N_VAR; ++j)
                dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];
            for (int j = 0; j < N_AUX; ++j)
                aux[storeix * nPts * N_AUX + j * nPts + i] = auxi[j];
        }
    }

    nStored[i] = storeix; //storeix ranged from 0 to nStored-1

    //write the final solution values to global memory.
    for (int j = 0; j < N_VAR; ++j)
        xf[j * nPts + i] = xi[j];

    // To get same RNG on repeat (non-continued) run, need to set the seed to same value
    for (int j = 0; j < N_RNGSTATE; ++j)
        RNGstate[j * nPts + i] = rd.state[j];

    // update dt to its final value (for adaptive stepper continue)
    d_dt[i] = dt;
}
