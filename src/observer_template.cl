// template to show function sigs

#ifndef OBSERVER_TEMPLATE_H
#define OBSERVER_TEMPLATE_H

//Number of features output. ObserverData must have a field F of this size, host also needs this info to allocate ObserverData and Feature arrays.
//TODO: could be defined by host as compiler flag have hardcoded value in only one place: host observer switch block.

#include "clODE_utilities.cl" //common math functions
#include "observers.cl"       //data structs

//typedef to select the correct observer struct from the header
typedef struct ObserverDataBasic
{
    int eventcount;
    int stepcount;
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

//prepare observer data to ensure it is ready for continuation if needed
//void prepareObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
//eg. times of last events need to be shifted left of t0. If last run had t0=0, then will be negative. However might continue from t0!=0, then need to add new t0.
//}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype *tspan)
{
    //eg. times of last events need to be shifted left of t0
}

#endif // OBSERVER_TEMPLATE_H
