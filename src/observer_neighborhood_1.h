//

#ifndef OBSERVER_NEIGHBORHOOD_1_H
#define OBSERVER_NEIGHBORHOOD_1_H

#include "observers.h"
#include "clODE_utilities.h"

//Use a solution buffer of 3 steps to detect extrema
#define BUFFER_SIZE 3

typedef struct ObserverData{
	realtype x0[N_VAR];
    int stepcount;
} ObserverData;


void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	od->stepcount=0;
}

//nothing to do - One pass detector
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {	
}

//all features are per-timestep
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[],  ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred) {
	++od->stepcount;
}

//event is return to neighborhood of initial/specified point in statespace
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
}

//no events
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	return false;
}

//no events
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	return false;
}

void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype * F, const int i, const int nPts) {
	int ix=0; 
	F[ix++*nPts+i]=od->xTrajectoryMax;
	F[ix++*nPts+i]=od->xTrajectoryMin;
	F[ix++*nPts+i]=od->xTrajectoryMean;
	F[ix++*nPts+i]=od->dxTrajectoryMax;
	F[ix++*nPts+i]=od->dxTrajectoryMin;
	F[ix++*nPts+i]=od->stepcount;
}


#endif // OBSERVER_NEIGHBORHOOD_1_H

