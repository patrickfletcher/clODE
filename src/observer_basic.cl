

#ifndef OBSERVER_BASIC_H
#define OBSERVER_BASIC_H

#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"

typedef struct ObserverData
{
    realtype xTrajectoryMax;
    realtype xTrajectoryMin;
    realtype xTrajectoryMean;
    realtype dxTrajectoryMax;
    realtype dxTrajectoryMin;
    int eventcount;
    int stepcount;
} ObserverData;

void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->xTrajectoryMax = -BIG_REAL;
    od->xTrajectoryMin = BIG_REAL;
    od->xTrajectoryMean = RCONST(0.0);
    od->dxTrajectoryMax = -BIG_REAL;
    od->dxTrajectoryMin = BIG_REAL;
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
    od->xTrajectoryMax = fmax(xi[op->fVarIx], od->xTrajectoryMax);
    od->xTrajectoryMin = fmin(xi[op->fVarIx], od->xTrajectoryMin);
    runningMean(&od->xTrajectoryMean, xi[op->fVarIx], od->stepcount);
    od->dxTrajectoryMax = fmax(dxi[op->fVarIx], od->dxTrajectoryMax);
    od->dxTrajectoryMin = fmin(dxi[op->fVarIx], od->dxTrajectoryMin);
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

void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype *F, const int i, const int nPts)
{
    int ix = 0;
    F[ix++ * nPts + i] = od->xTrajectoryMax;
    F[ix++ * nPts + i] = od->xTrajectoryMin;
    F[ix++ * nPts + i] = od->xTrajectoryMean;
    F[ix++ * nPts + i] = od->dxTrajectoryMax;
    F[ix++ * nPts + i] = od->dxTrajectoryMin;
    F[ix++ * nPts + i] = od->stepcount;
}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype tspan[])
{
    //nothing to do
}

#endif // OBSERVER_BASIC_H
