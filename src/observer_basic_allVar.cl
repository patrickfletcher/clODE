

#ifndef OBSERVER_BASIC_ALLVAR_H
#define OBSERVER_BASIC_ALLVAR_H

#include "clODE_utilities.cl"
#include "observers.cl"

typedef struct ObserverData
{
    realtype xTrajectoryMax[N_VAR];
    realtype xTrajectoryMin[N_VAR];
    realtype xTrajectoryMean[N_VAR];
    realtype dxTrajectoryMax[N_VAR];
    realtype dxTrajectoryMin[N_VAR];

    realtype auxTrajectoryMax[N_AUX];
    realtype auxTrajectoryMin[N_AUX];
    realtype auxTrajectoryMean[N_AUX];

    int eventcount;
    int stepcount;
} ObserverData;

void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->stepcount = 0;
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xTrajectoryMax[j] = -BIG_REAL;
        od->xTrajectoryMin[j] = BIG_REAL;
        od->xTrajectoryMean[j] = RCONST(0.0);
        od->dxTrajectoryMax[j] = -BIG_REAL;
        od->dxTrajectoryMin[j] = BIG_REAL;
    }
    for (int j = 0; j < N_AUX; ++j)
    {
        od->auxTrajectoryMax[j] = -BIG_REAL;
        od->auxTrajectoryMin[j] = BIG_REAL;
        od->auxTrajectoryMean[j] = RCONST(0.0);
    }
    od->eventcount = 0;
    od->stepcount = 0;
}

//nothing to do - One pass detector
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
}

//all features are per-timestep (eventOccurred unused)
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred)
{
    ++od->stepcount;
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xTrajectoryMax[j] = fmax(xi[j], od->xTrajectoryMax[j]);
        od->xTrajectoryMin[j] = fmin(xi[j], od->xTrajectoryMin[j]);
        runningMean(&od->xTrajectoryMean[j], xi[j], od->stepcount);
        od->dxTrajectoryMax[j] = fmax(dxi[j], od->dxTrajectoryMax[j]);
        od->dxTrajectoryMin[j] = fmin(dxi[j], od->dxTrajectoryMin[j]);
    }
    for (int j = 0; j < N_AUX; ++j)
    {
        od->auxTrajectoryMax[j] = fmax(auxi[j], od->auxTrajectoryMax[j]);
        od->auxTrajectoryMin[j] = fmin(auxi[j], od->auxTrajectoryMin[j]);
        runningMean(&od->auxTrajectoryMean[j], auxi[j], od->stepcount);
    }
}

//no events
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
}

//no events
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    return false;
}

//no events
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    return false;
}

void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype *F, int i, int nPts)
{
    //N_FEATS=5*N_VAR + 3*N_AUX + 1
    int ix = 0;
    for (int j = 0; j < N_VAR; ++j)
    {
        F[ix++ * nPts + i] = od->xTrajectoryMax[j]-od->xTrajectoryMin[j];
    }
    for (int j = 0; j < N_AUX; ++j)
    {
        F[ix++ * nPts + i] = od->auxTrajectoryMax[j]-od->auxTrajectoryMin[j];
    }
    
    for (int j = 0; j < N_VAR; ++j)
    {
        F[ix++ * nPts + i] = od->xTrajectoryMax[j];
        F[ix++ * nPts + i] = od->xTrajectoryMin[j];
        F[ix++ * nPts + i] = od->xTrajectoryMean[j];
        F[ix++ * nPts + i] = od->dxTrajectoryMax[j];
        F[ix++ * nPts + i] = od->dxTrajectoryMin[j];
    }
    for (int j = 0; j < N_AUX; ++j)
    {
        F[ix++ * nPts + i] = od->auxTrajectoryMax[j];
        F[ix++ * nPts + i] = od->auxTrajectoryMin[j];
        F[ix++ * nPts + i] = od->auxTrajectoryMean[j];
    }
    F[ix++ * nPts + i] = od->stepcount;
}
//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype *tspan)
{
    //nothing to do
}

#endif // OBSERVER_BASIC_H
