//

#ifndef OBSERVER_THRESHOLD_2_H
#define OBSERVER_THRESHOLD_2_H

#include "clODE_utilities.cl"
#include "observers.cl"

typedef struct ObserverData
{
    //compute using max/min operators
    realtype xTrajectoryMax;
    realtype xTrajectoryMin;
    realtype xTrajectoryAmp;

    realtype dxTrajectoryMax;
    realtype dxTrajectoryMin;
    realtype dxTrajectoryAmp;

    // realtype auxTrajectoryMax[N_AUX];
    // realtype auxTrajectoryMin[N_AUX];

    //thresholds
    realtype xUp;
    realtype xDown;
    realtype dxUp;
    realtype dxDown;

    //event data, update at event detection
    realtype tThisEvent;
    // realtype xThisEvent[N_VAR];
    // realtype dxThisEvent[N_VAR];
    // realtype auxThisEvent[N_AUX];
    //~ realtype nDistinctEvents[N_DISTINCT_MAX];

    //event data, update after computing event-triggered features
    realtype tLastEvent;
    // realtype xLastEvent[N_VAR];
    // realtype dxLastEvent[N_VAR];
    // realtype auxLastEvent[N_AUX];

    //compute using max/min operators, reset after computing event-triggered features
    // realtype xEventMax[N_VAR];
    // realtype xEventMin[N_VAR];
    // realtype dxEventMax[N_VAR];
    // realtype dxEventMin[N_VAR];
    // realtype auxEventMax[N_AUX];
    // realtype auxEventMin[N_AUX];

    realtype tAtDown;
    // realtype xAtDown[N_VAR];
    // realtype dxAtDown[N_VAR];
    // realtype auxAtDown[N_AUX];

    //local extrema (optionally computed, using sign change in slope)
    realtype tLastMax;
    // realtype xLastMax[N_VAR];
    // realtype auxLastMax[N_AUX];

    realtype tLastMin;
    // realtype xLastMin[N_VAR];
    // realtype auxLastMin[N_AUX];

    int stepcount;
    int eventcount;
    int mincount;
    int eventmincount;
    int maxcount;
    int eventmaxcount;
    bool inUpstate;
} ObserverData;

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
}

//restricted per-timestep update of observer data for initializing event detector
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->xTrajectoryMax[op->eVarIx] = fmax(xi[op->eVarIx], od->xTrajectoryMax[op->eVarIx]);
    od->xTrajectoryMin[op->eVarIx] = fmin(xi[op->eVarIx], od->xTrajectoryMin[op->eVarIx]);
    od->dxTrajectoryMax[op->eVarIx] = fmax(dxi[op->eVarIx], od->dxTrajectoryMax[op->eVarIx]);
    od->dxTrajectoryMin[op->eVarIx] = fmin(dxi[op->eVarIx], od->dxTrajectoryMin[op->eVarIx]);
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->xTrajectoryAmp = od->xTrajectoryMax - od->xTrajectoryMin;
    od->dxTrajectoryAmp = od->dxTrajectoryMax - od->dxTrajectoryMin;

    od->xUp = od->xTrajectoryMin + op->xUpThresh * od->xTrajectoryAmp;
	if(op->xDownThresh>RCONST(0.0))
		od->xDown = od->xTrajectoryMin + op->xDownThresh * od->xTrajectoryAmp;
	else
		obs->xDown = obs->xUp;

    od->dxUp = od->dxTrajectoryMin + op->dxUpThresh * od->dxTrajectoryAmp;
	if(op->xDownThresh>RCONST(0.0))
        od->dxDown = od->dxTrajectoryMin + op->dxDownThresh * od->dxTrajectoryAmp;
	else
		obs->dxDown = obs->dxUp;

	//determine if we are up or down.
	obs->inUpstate=(xi[op->eVarIx] > obs->xUp) ? true : false;
}

//per-timestep check for an event.  Option: refine event (t,x,dx,aux) within the timestep with interpolation
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    return false;
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred)
{
    if (eventOccurred)
    { //reset event-based quantities
    }
    else
    { //update per-step data and check for any intra-event special points of interest
    }
}

//Perform and post-integration cleanup and write desired features into the global array F
void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype *F, int i, int nPts)
{
    //Number of features is determined by this function. Must hardcode that number into the host program in order to allocate memory for F...
}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype tspan[])
{
    //nothing to do
}

#endif // OBSERVER_THRESHOLD_2_H
