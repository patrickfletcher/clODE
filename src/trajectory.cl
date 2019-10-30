//trajectory: stores in global variables directly

//TODO: optionally store only a subset of variables (allow bigger nPts): sp.varIx vector of indices? already in observerpars...???
//TODO: alternate storage at specified time points only - host sets t vector, interp and store x/dx/aux whenever ti passes t[nextstoreix]
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
    __global realtype *auxf,            //final value of aux variables 	[nPts*nAux]
    __global ulong *RNGstate,           //state for RNG				[nPts*nRNGstate]
    __global realtype *t,               //
    __global realtype *x,               //
    __global realtype *dx,              //
    __global realtype *aux,             //
    int nStoreMax,
    __global int *nStored)
{
    int i = get_global_id(0);
    int nPts = get_global_size(0);

    realtype ti, dt;
    realtype p[N_PAR], xi[N_VAR], dxi[N_VAR], auxi[N_AUX], wi[N_WIENER];
    rngData rd;
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
    for (int j = 0; j < N_WIENER; ++j)
        wi[j] = RCONST(0.0);
#else
    for (int j = 0; j < N_WIENER; ++j)
        wi[j] = randn(&rd) / sqrt(dt);
#endif

    getRHS(ti, xi, p, dxi, auxi, wi); //slope at initial point, needed for FSAL steppers (bs23, dorpri5) and for DX output

    //store the initial point

    int storeix = 0;
    t[storeix + i] = tspan[0];
    for (int j = 0; j < N_VAR; ++j)
        x[storeix * nPts * N_VAR + j * nPts + i] = xi[j];

    for (int j = 0; j < N_VAR; ++j)
        dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];

    for (int j = 0; j < N_AUX; ++j)
        aux[storeix * nPts * N_AUX + j * nPts + i] = auxi[j];

    //time-stepping loop, main time interval
    int step = 0;
    int stepflag = 0;
    while (ti < tspan[1] && step < sp->max_steps && stepflag == 0 && storeix < nStoreMax)
    {

        ++step;
#ifdef ADAPTIVE_STEPSIZE
        //leave the wi=0 for adaptive steppers
        stepflag = stepper(&ti, xi, dxi, p, sp, &dt, tspanPtr, auxi, wi);
#else
        //update Wiener variables - fixed size steppers can scale by dt here
        for (int j = 0; j < N_WIENER; ++j)
            wi[j] = randn(&rd) / sqrt(dt); //NOTE: divide by sqrt(dt) because Euler will multiply this by dt in the stepper.

        stepper(&ti, xi, dxi, p, dt, auxi, wi);
        ti = tspan[0] + step * dt;                                        //purify ti - Gets nSteps correct, but incompatible with shrinking final step without conditional to check if doing the last step
                                                                          //FSAL: dxi is at new ti, Not FSAL: dxi is at old ti
#endif

        //store every sp.nout'th step after the initial point
        if (step % sp->nout == 0)
        {
            ++storeix;

            t[storeix * nPts + i] = ti; //adaptive steppers give different timepoints for each trajectory

            for (int j = 0; j < N_VAR; ++j)
                x[storeix * nPts * N_VAR + j * nPts + i] = xi[j];

            for (int j = 0; j < N_VAR; ++j)
#ifdef FSAL_STEP_PROPERTY
                dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];
#else
                dx[(storeix - 1) * nPts * N_VAR + j * nPts + i] = dxi[j]; //without FSAL, dxi returned by stepper is evaluated at beginning of timestep
#endif

            for (int j = 0; j < N_AUX; ++j)
                aux[storeix * nPts * N_AUX + j * nPts + i] = auxi[j];
        }
    }
    //~ aux[i]=0;

#if !defined(FSAL_STEP_PROPERTY)
    //get the slope at the last point
    #ifdef ADAPTIVE_STEPSIZE
        // for (int j = 0; j < N_WIENER; ++j)
        //     wi[j] = RCONST(0.0);
    #else
        for (int j = 0; j < N_WIENER; ++j)
            wi[j] = randn(&rd) / sqrt(dt);
    #endif

    getRHS(ti, xi, p, dxi, auxi, wi);

    for (int j = 0; j < N_VAR; ++j)
        dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];
#endif

    //get device arrays ready to continue
    for (int j = 0; j < N_VAR; ++j)
        xf[j * nPts + i] = xi[j];

    for (int j = 0; j < N_AUX; ++j)
        auxf[j * nPts + i] = auxi[j];

    for (int j = 0; j < N_RNGSTATE; ++j)
        RNGstate[j * nPts + i] = rd.state[j];

    nStored[i] = storeix; //storeix ranged from 0 to nStored-1
}
