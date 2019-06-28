//

#ifndef OBSERVER_NEIGHBORHOOD_2_H
#define OBSERVER_NEIGHBORHOOD_2_H

#include "observers.h"
#include "clODE_utilities.h"

#define TWO_PASS_EVENT_DETECTOR

//typedef to select the correct observer struct from the header
typedef struct ObserverData{
	int eventcount;
    int stepcount;
    int buffer_filled;

	realtype tbuffer[3];
    realtype xbuffer[3]; 
    realtype dxbuffer[3];

	realtype x0[N_VAR]; //point for neighborhood return
	realtype xLast[N_VAR];
	realtype xMaxDelta[N_VAR]; //store the maximum delta-x = xi[j]-xLast[j] per time step in each dimension (for size of epsilon "ball")

	realtype xTrajectoryMax[N_VAR]; 
	realtype xTrajectoryMin[N_VAR];
	realtype xRange[N_VAR];
	realtype dxTrajectoryMax[N_VAR];
	realtype dxTrajectoryMin[N_VAR];
	realtype dxRange[N_VAR];

	bool isInNhood;


} ObserverData;

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {	
	
	od->eventcount=0;
	od->stepcount=0;
	od->buffer_filled=0;
	od->isInNhood=false;

	for (int j=0; j<N_VAR;++j) {
		od->xLast[j]=xi[j];
		od->xMaxDelta[j]=-BIG_REAL;

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
		od->xMaxDelta[j]=MAX(xi[j]-od->xLast[j], od->xMaxDelta);
		od->xLast[j]=xi[j];

		od->xTrajectoryMax[j]=MAX(xi[j], od->xTrajectoryMax[j]);
		od->xTrajectoryMin[j]=MIN(xi[j], od->xTrajectoryMin[j]);
		od->dxTrajectoryMax[j]=MAX(dxi[j], od->dxTrajectoryMax[j]);
		od->dxTrajectoryMin[j]=MIN(dxi[j], od->dxTrajectoryMin[j]);
	}
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	
	//find slowest variable (min dx amplitude?)

    for (int j=0; j<N_VAR; ++j) {
		od->x0[j]=xi[j];
        od->xRange[j]=od->xTrajectoryMax[j]-od->xTrajectoryMin[j];
        od->dxRange[j]=od->dxTrajectoryMax[j]-od->dxTrajectoryMin[j];
    }
	od->isInNhood = true; //we are at the center of the neighborhood, x0.
	
}

//check for entry into epsilon ball surrounding od->x0
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	if (od->buffer_filled==1) {
		bool lastWasInNhood = od->isInNhood;
    	for (int j=0; j<N_VAR; ++j) {
			od->isInNhood=od->isInNhood & fabs(od->x0[j]-xi[j])<2*od->xMaxDelta[j];
		}
		return (od->isInNhood & ~lastWasInNhood)
		}
	}
	return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
// - get (t,x) at local max, compute: IMI, tMaxMin, amp
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	realtype tThisMax, xThisMax;
	realtype thisIMI, thisTMaxMin, thisAmp;
	
	++od->eventcount;
	
	
	return false; //not terminal 
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
// - advance solution/slope buffers
// - check for any intermediate special points & store their info
// - reset intermediates upon local max event detection
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[],  ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred) {
	
    ++od->stepcount;
    
    if (od->stepcount>3 && od->buffer_filled==0) { od->buffer_filled=1;}
    
	//advance solution buffer
	od->tbuffer[0]=od->tbuffer[1];
	od->tbuffer[1]=od->tbuffer[2];
	od->tbuffer[2]=*ti;	

	for (int i=0; i<N_VAR; ++i) {
		od->xbuffer[i*3+0]=od->xbuffer[i*3+1];
		od->xbuffer[i*3+1]=od->xbuffer[i*3+2];
		od->xbuffer[i*3+2]=yi[i];	
	}
	
	od->dxbuffer[0]=od->dxbuffer[1];
	od->dxbuffer[1]=od->dxbuffer[2];	
	od->dxbuffer[2]=dyi[op->varIx];
	
	//global dxMax, dxMin
	od->dxMax=MAX(od->dxMax,dxi[op->fVarIx]); 
	od->dxMin=MIN(od->dxMin,dxi[op->fVarIx]);
	runningMean(&od->xTrajectoryMean, xi[op->fVarIx], od->stepcount);
	
	
	//reset intermediate feature storage. Should this be here, or in "computeEventFeatures"?
	if (eventOccurred) { 
		
	}
}

//Perform and post-integration cleanup and write desired features into the global array F
void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype * F, int i, int nPts) {
	int ix=0;
	F[ix++*nPts+i]=od->eventcount>1?od->IMI[0]:0;
	F[ix++*nPts+i]=od->eventcount>1?od->IMI[1]:0;
	F[ix++*nPts+i]=od->IMI[2];
	F[ix++*nPts+i]=od->eventcount>1?od->amp[0]:0;
	F[ix++*nPts+i]=od->eventcount>1?od->amp[1]:0;
	F[ix++*nPts+i]=od->amp[2];
	//F[ix++*nPts+i]=od->eventcount>1?od->tMaxMin[0]:0;
	//F[ix++*nPts+i]=od->eventcount>1?od->tMaxMin[1]:0;
	//F[ix++*nPts+i]=od->tMaxMin[2];
	F[ix++*nPts+i]=od->eventcount>0?od->xMax[0]:xi[op->fVarIx];
	F[ix++*nPts+i]=od->eventcount>0?od->xMax[1]:xi[op->fVarIx];
	F[ix++*nPts+i]=od->eventcount>0?od->xMax[2]:xi[op->fVarIx];
	F[ix++*nPts+i]=od->eventcount>0?od->xMin[0]:xi[op->fVarIx];
	F[ix++*nPts+i]=od->eventcount>0?od->xMin[1]:xi[op->fVarIx];
	F[ix++*nPts+i]=od->eventcount>0?od->xMin[2]:xi[op->fVarIx];
	F[ix++*nPts+i]=od->dxMax;
	F[ix++*nPts+i]=od->dxMin;
	F[ix++*nPts+i]=od->xTrajectoryMean;
	F[ix++*nPts+i]=od->eventcount;
	F[ix++*nPts+i]=od->stepcount;
}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype tspan[]) {
	//shift all time-based observer members left by [tf-t0]
	realtype T=*ti-tspan[0];
	od->tLastMax=od->tLastMax-T;
	od->tLastMin=od->tLastMin-T;
	for (int i=0; i<BUFFER_SIZE; ++i) {
		od->tbuffer[i]=od->tbuffer[i]-T; }
}


#endif // OBSERVER_NEIGHBORHOOD_2_H

