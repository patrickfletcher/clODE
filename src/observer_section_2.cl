//

#ifndef OBSERVER_SECTION_2_H
#define OBSERVER_SECTION_2_H

#include "clODE_utilities.cl"
#include "observers.cl"

//Use a solution buffer of 3 steps to detect extrema
#define BUFFER_SIZE 3

typedef struct ObserverData
{

    //compute using max/min operators
    realtype xTrajectoryMax[N_VAR];
    realtype xTrajectoryMin[N_VAR];
    realtype dxTrajectoryMax[N_VAR];
    realtype dxTrajectoryMin[N_VAR];
    realtype auxTrajectoryMax[N_AUX];
    realtype auxTrajectoryMin[N_AUX];

    //event data, update at event detection
    realtype tThisEvent;
    realtype xThisEvent[N_VAR];
    realtype dxThisEvent[N_VAR];
    realtype auxThisEvent[N_AUX];
    //~ realtype nDistinctEvents[N_DISTINCT_MAX];

    //event data, update after computing event-triggered  features
    realtype tLastEvent;
    realtype xLastEvent[N_VAR];
    realtype dxLastEvent[N_VAR];
    realtype auxLastEvent[N_AUX];

    //compute using max/min operators, reset after computing event-triggered features
    realtype xEventMax[N_VAR];
    realtype xEventMin[N_VAR];
    realtype dxEventMax[N_VAR];
    realtype dxEventMin[N_VAR];
    realtype auxEventMax[N_AUX];
    realtype auxEventMin[N_AUX];

    //thresholds
    realtype xUp;
    realtype xDown;
    realtype dxUp;
    realtype dxDown;

    realtype tAtDown;
    realtype xAtDown[N_VAR];
    realtype dxAtDown[N_VAR];
    realtype auxAtDown[N_AUX];

    //local extrema (optionally computed, using sign change in slope)
    realtype tLastMax;
    realtype xLastMax[N_VAR];
    realtype dxLastMax[N_VAR];
    realtype auxLastMax[N_AUX];

    realtype tLastMin;
    realtype xLastMin[N_VAR];
    realtype dxLastMin[N_VAR];
    realtype auxLastMin[N_AUX];

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
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
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

#endif // OBSERVER_SECTION_2_H
