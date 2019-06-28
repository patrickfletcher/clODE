//use local maxima as events

#ifndef OBSERVER_LOCAL_EXTREMUM_H
#define OBSERVER_LOCAL_EXTREMUM_H

#include "observers.h"
#include "clODE_utilities.h"


//events are triggered at local maxima in the variable specified by op.fVarIx

//Use a solution buffer of 3 steps to detect extrema
#define BUFFER_SIZE 3

typedef struct ObserverData {
	int eventcount;
    int stepcount;
    int buffer_filled;

	realtype tbuffer[BUFFER_SIZE];
    realtype xbuffer[BUFFER_SIZE]; 
    realtype dxbuffer[BUFFER_SIZE];
    
	realtype xTrajectoryMean;
	//~ realtype nDistinctEvents[N_DISTINCT_MAX];
	
	//local max data
	realtype tLastMax;
	
	//one minimum between two maxima
	realtype tLastMin;
	realtype xLastMin;
	
	//record extrema in dx and aux too
	realtype dxMax;
	realtype dxMin;
	
	realtype IMI[3]; //max/min/mean
	realtype amp[3]; //max/min/mean
	//realtype tMaxMin[3]; //max/min/mean
	realtype xMax[3]; //max/min/mean
	realtype xMin[3]; //max/min/mean

}ObserverData;
//size: (3*BUFFER_SIZE + 6 + 3*3)*realtype + 2 int 


//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	
	od->dxMax=-BIG_REAL;
	od->dxMin= BIG_REAL;
	
	od->IMI[0]=-BIG_REAL;
	od->IMI[1]= BIG_REAL;
	od->IMI[2]= 0;
	
	od->amp[0]=-BIG_REAL;
	od->amp[1]= BIG_REAL;
	od->amp[2]= 0;
	
	//od->tMaxMin[0]=-BIG_REAL;
	//od->tMaxMin[1]= BIG_REAL;
	//od->tMaxMin[2]= 0;
	
	
	od->xMax[0]=-BIG_REAL;
	od->xMax[1]= BIG_REAL;
	od->xMax[2]= 0;
	
	od->xMin[0]=-BIG_REAL;
	od->xMin[1]= BIG_REAL;
	od->xMin[2]= 0;
	
	od->xLastMin=BIG_REAL;
	od->tLastMin=0;

	od->eventcount=0;
	od->stepcount=0;
	od->buffer_filled=0;
}

//no warmup needed
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	//nothing to do
}

//check buffer of slopes for local max
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	if (od->buffer_filled==1) {
		return (od->dxbuffer[1] >= -op->eps_dx && od->dxbuffer[2] < -op->eps_dx );
	}
	return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
// - get (t,x) at local max, compute: IMI, tMaxMin, amp
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op) {
	realtype tThisMax, xThisMax;
	realtype thisIMI, thisTMaxMin, thisAmp;
	
	++od->eventcount;
	
	//Simple Max of xbuffer  
	int ix=0;
	maxOfArray(od->xbuffer, BUFFER_SIZE, &xThisMax, &ix);
	tThisMax=od->tbuffer[ix];
	
	//Quadratic interpolation for improved accuracy BROKEN??
	//quadraticInterpVertex(od->tbuffer, od->xbuffer, &tThisMax, &xThisMax);
	
	od->xMax[0]=MAX(xThisMax, od->xMax[0]);
	od->xMax[1]=MIN(xThisMax, od->xMax[1]);
	runningMean(&od->xMax[2],xThisMax,od->eventcount-1);
	
	od->xMin[0]=MAX(od->xLastMin, od->xMin[0]);
	od->xMin[1]=MIN(od->xLastMin, od->xMin[1]);
	runningMean(&od->xMin[2],od->xLastMin,od->eventcount-1);
	
	if (od->eventcount>1) {//implies tLastMax, tLastMin and xLastMin are set
		
		//max/min/mean IMI
		thisIMI = tThisMax - od->tLastMax;
		//if(thisIMI > op->minIMI) {
			od->IMI[0]=MAX(thisIMI, od->IMI[0]);
			od->IMI[1]=MIN(thisIMI, od->IMI[1]);
			runningMean(&od->IMI[2],thisIMI,od->eventcount-1);
		//}
		
		//max/min/mean tMaxMin
		//thisTMaxMin=tThisMax-od->tLastMin;
		//od->tMaxMin[0]=MAX(thisTMaxMin, od->tMaxMin[0]);
		//od->tMaxMin[1]=MIN(thisTMaxMin, od->tMaxMin[1]);
		//runningMean(&od->tMaxMin[2],thisTMaxMin,od->eventcount-1);
		
		//max/min/mean amp
		thisAmp=xThisMax-od->xLastMin;
		od->amp[0]=MAX(thisAmp, od->amp[0]);
		od->amp[1]=MIN(thisAmp, od->amp[1]);
		runningMean(&od->amp[2],thisAmp,od->eventcount-1);
	}
	
	//update stored "last" values
	od->tLastMax=tThisMax;
	
	return false; //not terminal 
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
// - advance solution/slope buffers
// - check for any intermediate special points & store their info
// - reset intermediates upon local max event detection
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[],  ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred) {
	
    ++od->stepcount;
    
    if (od->stepcount>BUFFER_SIZE && od->buffer_filled==0) { od->buffer_filled=1;}
    
	//advance solution buffer
	for (int i=0; i<BUFFER_SIZE-1; ++i) {
		od->tbuffer[i]=od->tbuffer[i+1];
		od->xbuffer[i]=od->xbuffer[i+1];
		od->dxbuffer[i]=od->dxbuffer[i+1];
	}
	od->tbuffer[BUFFER_SIZE-1]=*ti;
	od->xbuffer[BUFFER_SIZE-1]=xi[op->fVarIx];
	od->dxbuffer[BUFFER_SIZE-1]=dxi[op->fVarIx];
	
	//global dxMax, dxMin
	od->dxMax=MAX(od->dxMax,dxi[op->fVarIx]); 
	od->dxMin=MIN(od->dxMin,dxi[op->fVarIx]);
	runningMean(&od->xTrajectoryMean, xi[op->fVarIx], od->stepcount);
	
	//local min check - one between each max - simply overwrite tLastMin, xLastMin
	// if (od->dxbuffer[1] <= op->eps_dx && od->dxbuffer[2] > op->eps_dx ) {
	// 	int ix=0;
	// 	minOfArray(od->xbuffer, BUFFER_SIZE, &od->xLastMin, &ix);
	// 	od->tLastMin=od->tbuffer[ix];
	// }
	od->xLastMin=MIN(od->xLastMin,xi[op->fVarIx]);
	
	//reset intermediate feature storage. Should this be here, or in "computeEventFeatures"?
	if (eventOccurred) { 
		//od->dxMax=-BIG_REAL; //only if looking for event dxMax/min
		//od->dxMin= BIG_REAL;
		od->xLastMin=BIG_REAL;
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

#endif //OBSERVER_LOCAL_EXTREMUM_H
