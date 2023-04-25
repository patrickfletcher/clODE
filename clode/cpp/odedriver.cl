#include "clODE_random.cl"
#include "clODE_struct_defs.cl"
#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"
#include "steppers.cl"

__kernel void odedriver(
    __constant realtype *tspan,         //time vector [t0,tf] - adds (tf-t0) to these at the end
    __global realtype *x0,              //initial state 				[nPts*nVar]
    __constant realtype *pars,          //parameter values				[nPts*nPar]
    __constant struct SolverParams *sp, //dtmin/max, tols, etc
    __global realtype *xf,              //final state 				[nPts*nVar]
    __global ulong *RNGstate,            //state for RNG					[nPts*nRNGstate]
    __global ObserverData *OData,        //Observer data. Assume it is initialized externally by initialize observer kernel!
    __constant struct ObserverParams *opars,
    __global realtype *F,  //feature results
    __global realtype *t,  // trajectory storage: t, x, dx, aux
    __global realtype *x,
    __global realtype *dx,
    __global realtype *aux,
    __global int *nStored)
{

    int i = get_global_id(0);
    int nPts = get_global_size(0);

    realtype ti, dt;
    realtype p[N_PAR], xi[N_VAR], dxi[N_VAR], auxi[N_AUX], wi[N_WIENER];
    rngData rd;
    __constant realtype * const tspanPtr = tspan;

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
	getRHS(ti, xi, p, dxi, auxi, wi); //slope at initial point, needed for FSAL steppers (bs23, dorpri5)

    //store the initial point --> could be set elsewhere using an "init" kernel?
    int storeix = 0;
    t[storeix + i] = tspan[0];
    for (int j = 0; j < N_VAR; ++j)
        x[storeix * nPts * N_VAR + j * nPts + i] = xi[j];

    for (int j = 0; j < N_VAR; ++j)
        dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];

    for (int j = 0; j < N_AUX; ++j)
        aux[storeix * nPts * N_AUX + j * nPts + i] = auxi[j];

    ObserverData odata = OData[i]; //private copy of observer data. Assume it is initialized externally by initialize observer kernel

    //time-stepping loop, main time interval
    int step = 0;
    int stepflag = 0;
    bool eventOccurred;
    bool terminalEvent;
    while (ti < tspan[1] && step < sp->max_steps && storeix < sp->max_store)
    {

        ++step;
        ++odata.stepcount;
        stepflag = stepper(&ti, xi, dxi, p, sp, &dt, tspanPtr, auxi, wi, step, &rd);
        if (stepflag!=0)
            break;

        eventOccurred = eventFunction(&ti, xi, dxi, auxi, &odata, opars);
        if (eventOccurred)
        {
            terminalEvent = computeEventFeatures(&ti, xi, dxi, auxi, &odata, opars);
            if (terminalEvent)
            {
                break;
            };
        }

        updateObserverData(&ti, xi, dxi, auxi, &odata, opars); //TODO: if not FSAL, dxi buffer is delayed by one. (dxi is slope at LAST timestep)

        //store every sp.nout'th step after the initial point
        if (step % sp->nout == 0)
        {
            ++storeix;

            t[storeix * nPts + i] = ti; //adaptive steppers give different timepoints for each trajectory

            for (int j = 0; j < N_VAR; ++j)
                x[storeix * nPts * N_VAR + j * nPts + i] = xi[j];

            for (int j = 0; j < N_VAR; ++j)
                dx[storeix * nPts * N_VAR + j * nPts + i] = dxi[j];

            for (int j = 0; j < N_AUX; ++j)
                aux[storeix * nPts * N_AUX + j * nPts + i] = auxi[j];
        }
    }

    //readout features of interest and write to global F:
    finalizeFeatures(&ti, xi, dxi, auxi, &odata, opars, F, i, nPts);

    //finalize observerdata for possible continuation
    finalizeObserverData(&ti, xi, dxi, auxi, &odata, opars, tspan);


    //write the final solution values to global memory.
    for (int j = 0; j < N_VAR; ++j)
        xf[j * nPts + i] = xi[j];

    // To get same RNG on repeat (non-continued) run, need to set the seed to same value
    for (int j = 0; j < N_RNGSTATE; ++j)
        RNGstate[j * nPts + i] = rd.state[j];

    OData[i] = odata;
    
    nStored[i] = storeix; //storeix ranged from 0 to nStored-1
}
