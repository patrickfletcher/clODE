//

#ifndef OBSERVER_NEIGHBORHOOD_2_H
#define OBSERVER_NEIGHBORHOOD_2_H

#include "observers.h"
#include "clODE_utilities.h"

//typedef to select the correct observer struct from the header
typedef struct ObserverData{
	realtype x0[N_VAR]; //point for neighborhood return
	realtype xTrajectoryMax[N_VAR]; 
	realtype xTrajectoryMin[N_VAR];
	realtype dxTrajectoryMax[N_VAR];
	realtype dxTrajectoryMin[N_VAR];
    int stepcount;
} ObserverData;

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {	od->stepcount=0;
	for (int j=0; j<N_VAR;++j) {
		od->xTrajectoryMax[j]=-BIG_REAL;
		od->xTrajectoryMin[j]= BIG_REAL;
		od->dxTrajectoryMax[j]=-BIG_REAL;
		od->dxTrajectoryMin[j]= BIG_REAL;
	}
}

//restricted per-timestep update of observer data for initializing event detector
// - get extent of trajectory in state space, and max/min slopes
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	for (int j=0; j<N_VAR;++j) {
		od->xTrajectoryMax[j]=MAX(xi[j],od->xTrajectoryMax[j]);
		od->xTrajectoryMin[j]=MIN(xi[j],od->xTrajectoryMin[j]);
		od->dxTrajectoryMax[j]=MAX(dxi[j],od->dxTrajectoryMax[j]);
		od->dxTrajectoryMin[j]=MIN(dxi[j],od->dxTrajectoryMin[j]);
	}
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	
	//find slowest variable (min dx amplitude)
	
}

//per-timestep check for an event.  Option: refine event (t,x,dx,aux) within the timestep with interpolation
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	return false; 
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[],  ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred) {
	if (eventOccurred) { //reset event-based quantities
	}
	else { //update per-step data and check for any intra-event special points of interest
	}
}

//Perform and post-integration cleanup and write desired features into the global array F
void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype * F, int i, int nPts) {
	//Number of features is determined by this function. Must hardcode that number into the host program in order to allocate memory for F... 
}


//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype tspan[]) {
	//nothing to do
}


#endif // OBSERVER_NEIGHBORHOOD_2_H

