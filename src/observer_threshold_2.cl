//

#ifndef OBSERVER_THRESHOLD_2_H
#define OBSERVER_THRESHOLD_2_H

#include "clODE_utilities.cl"
#include "observers.cl"

typedef struct ObserverData
{
    realtype tbuffer[3];
    realtype xbuffer[3];
    realtype dxbuffer[3];

    realtype xTrajectoryMean;
    realtype xGlobalMax;
    realtype xGlobalMin;
    realtype dxGlobalMax;
    realtype dxGlobalMin;

    //thresholds
    realtype xUp;
    realtype xDown;
    realtype dxUp;
    realtype dxDown;

    realtype tLastUp;
    realtype tThisDown;

    realtype tLastMax; //compare to tThis
    realtype tLastMin;
    realtype xLastMax;
    realtype xLastMin; 

    realtype nMaxima[3]; //max/min/mean
    realtype period[3];  //max/min/mean
    realtype upDuration[3];  //max/min/mean
    realtype downDuration[3];  //max/min/mean

    int thisNMaxima;
    int stepcount;
    int eventcount;
    bool buffer_filled;
    bool inUpstate;
} ObserverData;

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
}

//restricted per-timestep update of observer data for initializing event detector
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->xGlobalMax = fmax(od->xGlobalMax, xi[op->fVarIx]);
    od->xGlobalMin = fmin(od->xGlobalMin, xi[op->fVarIx]);
    od->dxGlobalMax = fmax(od->dxGlobalMax, dxi[op->fVarIx]);
    od->dxGlobalMin = fmin(od->dxGlobalMin, dxi[op->fVarIx]);
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    realtype xTrajectoryAmp = od->xGlobalMax - od->xGlobalMin;
    realtype dxTrajectoryAmp = od->dxGlobalMax - od->dxGlobalMin;

    od->xUp = od->xGlobalMin + op->xUpThresh * xTrajectoryAmp;
	if(op->xDownThresh>RCONST(0.0))
		od->xDown = od->xGlobalMin + op->xDownThresh * xTrajectoryAmp;
	else
		obs->xDown = obs->xUp;

    od->dxUp = od->dxGlobalMin + op->dxUpThresh * dxTrajectoryAmp;
	if(op->xDownThresh>RCONST(0.0))
        od->dxDown = od->dxGlobalMin + op->dxDownThresh * dxTrajectoryAmp;
	else
		obs->dxDown = obs->dxUp;

	//determine if we are up or down.
	obs->inUpstate=(xi[op->eVarIx] > obs->xUp) ? true : false;
}

//per-timestep check for an event.  Option: refine event (t,x,dx,aux) within the timestep with interpolation
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    //event is marked by upward threshold crossing
    if (xi[op->eVarIx] >= obs->xUp && dxi[op->eVarIx] >= obs->dxUp && !obs->inUpstate)
    {
        obs->inUpstate = true;
    }

    return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    realtype tThisUp;

    ++od->eventcount;

    //simplest: time and state of trajectory point that landed in the neighborhood.
    // tThisUp=*ti;

    //alts: linearly interp to get more accurate tThisUp within the last time step: xi, xi-1, ti, ti-1
    tThisUp = linearInterp(od->xbuffer[2], xi[op->eVarIx], od->tbuffer[2], *ti, op->xUp);

    if (od->eventcount > 1)
    {
        od->nMaxima[0] = fmax((realtype)od->thisNMaxima, od->nMaxima[0]); //cast to realtype (for mean)
        od->nMaxima[1] = fmin((realtype)od->thisNMaxima, od->nMaxima[1]);
        runningMean(&od->nMaxima[2], (realtype)od->thisNMaxima, od->eventcount - 1);

        realtype thisPeriod = tThisUp - od->tLastUp;
        od->period[0] = fmax(thisPeriod, od->period[0]);
        od->period[1] = fmin(thisPeriod, od->period[1]);
        runningMean(&od->period[2], thisPeriod, od->eventcount - 1);
        
        realtype thisUpDuration = od->tThisDown - od->tLastUp;
        od->upDuration[0] = fmax(thisUpDuration, od->upDuration[0]);
        od->upDuration[1] = fmin(thisUpDuration, od->upDuration[1]);
        runningMean(&od->upDuration[2], thisUpDuration, od->eventcount - 1);
        
        realtype thisDownDuration = tThisUp - od->tThisDown;
        od->period[0] = fmax(thisPeriod, od->period[0]);
        od->period[1] = fmin(thisPeriod, od->period[1]);
        runningMean(&od->period[2], thisPeriod, od->eventcount - 1);
    }

    od->tLastUp = tThisUp;
    od->thisNMaxima = 0;

    return false;
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred)
{
    //global
    od->xGlobalMax = fmax(od->xGlobalMax, xi[op->fVarIx]);
    od->xGlobalMin = fmin(od->xGlobalMin, xi[op->fVarIx]);
    od->dxGlobalMax = fmax(od->dxGlobalMax, dxi[op->fVarIx]);
    od->dxGlobalMin = fmin(od->dxGlobalMin, dxi[op->fVarIx]);
    runningMean(&od->xTrajectoryMean, xi[op->fVarIx], od->stepcount);
    if (od->xGlobalMax - od->xGlobalMin < op->minXamp)
    {
        return;
    }

    //advance solution buffer
    od->tbuffer[0] = od->tbuffer[1];
    od->tbuffer[1] = od->tbuffer[2];
    od->tbuffer[2] = *ti;
    od->xbuffer[0] = od->xbuffer[1];
    od->xbuffer[1] = od->xbuffer[2];
    od->xbuffer[2] = xi[op->fVarIx];
    od->dxbuffer[0] = od->dxbuffer[1];
    od->dxbuffer[1] = od->dxbuffer[2];
    od->dxbuffer[2] = dxi[op->fVarIx];

    if (od->stepcount > 2 && od->buffer_filled == 0)
    {
        od->buffer_filled = 1;
    }

    //buffer-based computations
    if (od->buffer_filled == 1)
    {
        
        //local max check in fVarIx - one between each min - simply overwrite tLastMax, xLastMax
        if (od->dxbuffer[1] > -op->eps_dx && od->dxbuffer[2] < -op->eps_dx)
        {
            od->thisNMaxima++;
        }

        //Check for downward threshold crossing (end of "active phase")
        if (xi[op->eVarIx] <= obs->xDown && dxi[op->eVarIx] <= obs->dxDown && obs->inUpstate)
        {
            obs->inUpstate = false;
        }
    }


    if (eventOccurred)
    { //reset event-based quantities
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
