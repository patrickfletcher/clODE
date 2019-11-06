//

#ifndef OBSERVER_THRESHOLD_2_H
#define OBSERVER_THRESHOLD_2_H

#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"

#define TWO_PASS_EVENT_DETECTOR

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

    realtype tLastEvent;
    realtype tThisDown;

    // realtype tLastMax; //compare to tThis
    // realtype tLastMin;
    // realtype xLastMax;
    // realtype xLastMin; 

    realtype nMaxima[3]; //max/min/mean
    realtype period[3];  //max/min/mean
    realtype upDuration[3];  //max/min/mean
    realtype downDuration[3];  //max/min/mean

    realtype maxDt;
    realtype minDt;

    int thisNMaxima;
    int stepcount;
    int eventcount;
    bool buffer_filled;
    bool inUpstate;
} ObserverData;

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->tbuffer[2] = *ti;
    od->xbuffer[2] = xi[op->fVarIx];
    od->dxbuffer[2] = dxi[op->fVarIx];

    od->xTrajectoryMean = RCONST(0.0);
    od->xGlobalMax = -BIG_REAL;
    od->xGlobalMin = BIG_REAL;
    od->dxGlobalMax = -BIG_REAL;
    od->dxGlobalMin = BIG_REAL;

    // od->xUp = RCONST(0.0);
    // od->xDown = RCONST(0.0);
    // od->dxUp = RCONST(0.0);
    // od->dxDown = RCONST(0.0);
    // od->tLastEvent = RCONST(0.0);
    // od->tThisDown = RCONST(0.0);
    // od->tLastMax = RCONST(0.0);
    // od->tLastMin = RCONST(0.0);
    // od->xLastMax = RCONST(0.0);
    // od->xLastMin = RCONST(0.0); 

    od->nMaxima[0] = -BIG_REAL;
    od->nMaxima[1] = BIG_REAL;
    od->nMaxima[2] = RCONST(0.0);

    od->period[0] = -BIG_REAL;
    od->period[1] = BIG_REAL;
    od->period[2] = RCONST(0.0);

    od->upDuration[0] = -BIG_REAL;
    od->upDuration[1] = BIG_REAL;
    od->upDuration[2] = RCONST(0.0);

    od->downDuration[0] = -BIG_REAL;
    od->downDuration[1] = BIG_REAL;
    od->downDuration[2] = RCONST(0.0);

    od->maxDt = -BIG_REAL;
    od->minDt = BIG_REAL;

    od->thisNMaxima=0;
    od->stepcount=0;
    od->eventcount=0;
    od->buffer_filled=false;
    od->inUpstate=false;
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
    //threshold in x
    realtype xTrajectoryAmp = od->xGlobalMax - od->xGlobalMin;
    od->xUp = od->xGlobalMin + op->xUpThresh * xTrajectoryAmp;
	if(op->xDownThresh>RCONST(0.0))
		od->xDown = od->xGlobalMin + op->xDownThresh * xTrajectoryAmp;
	else
		od->xDown = od->xUp;

    //threshold for dx
    od->dxUp =  op->dxUpThresh * od->dxGlobalMax;
	if(op->dxDownThresh>RCONST(0.0))
        od->dxDown = op->dxDownThresh * od->dxGlobalMin;
	else
		od->dxDown = od->dxGlobalMin;

	//determine if we are up or down.
	od->inUpstate=(xi[op->fVarIx] > od->xUp) ? true : false;
    // od->xGlobalMax = -BIG_REAL; //TODO: need to keep a separate "globalAmp" from init, or don't wipe globalMax/Min: for minAmp to work
    // od->xGlobalMin = BIG_REAL;
    // od->dxGlobalMax = -BIG_REAL;
    // od->dxGlobalMin = BIG_REAL;
}

//per-timestep check for an event.  Option: refine event (t,x,dx,aux) within the timestep with interpolation
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    //event is marked by upward threshold crossing
    return (xi[op->fVarIx] >= od->xUp && dxi[op->fVarIx] >= od->dxUp && !od->inUpstate);
    
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    realtype tThisUp;

    ++od->eventcount;
    od->inUpstate = true;

    //simplest: time and state of trajectory point that landed in the neighborhood.
    // tThisUp=*ti;

    //alts: linearly interp to get more accurate tThisUp within the last time step: xi-1, xi, ti-1, ti; get t @ x=xUp
    tThisUp = linearInterp(od->xbuffer[2], xi[op->fVarIx], od->tbuffer[2], *ti, od->xUp);

    if (od->eventcount > 1)
    {
        od->nMaxima[0] = fmax((realtype)od->thisNMaxima, od->nMaxima[0]); //cast to realtype (for mean)
        od->nMaxima[1] = fmin((realtype)od->thisNMaxima, od->nMaxima[1]);
        runningMean(&od->nMaxima[2], (realtype)od->thisNMaxima, od->eventcount - 1);

        realtype thisPeriod = tThisUp - od->tLastEvent;
        od->period[0] = fmax(thisPeriod, od->period[0]);
        od->period[1] = fmin(thisPeriod, od->period[1]);
        runningMean(&od->period[2], thisPeriod, od->eventcount - 1);
        
        realtype thisUpDuration = od->tThisDown - od->tLastEvent;
        od->upDuration[0] = fmax(thisUpDuration, od->upDuration[0]);
        od->upDuration[1] = fmin(thisUpDuration, od->upDuration[1]);
        runningMean(&od->upDuration[2], thisUpDuration, od->eventcount - 1);
        
        realtype thisDownDuration = tThisUp - od->tThisDown;
        od->downDuration[0] = fmax(thisDownDuration, od->downDuration[0]);
        od->downDuration[1] = fmin(thisDownDuration, od->downDuration[1]);
        runningMean(&od->downDuration[2], thisDownDuration, od->eventcount - 1);
    }

    od->tLastEvent = tThisUp;
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

    // if (od->xGlobalMax - od->xGlobalMin < op->minXamp)
    // {
    //     return;
    // }

    if (od->xGlobalMax - od->xGlobalMin > op->minXamp)
    {
        if (od->stepcount > 2 && od->buffer_filled == 0)
        {
            od->buffer_filled = 1;
        }

        //buffer-based computations
        if (od->buffer_filled == 1)
        {        
            //local max check in fVarIx - one between each min - simply overwrite tLastMax, xLastMax
            if (od->dxbuffer[1] >= -op->eps_dx && od->dxbuffer[2] < -op->eps_dx)
            {
                od->thisNMaxima++;
                // int ix;
                // maxOfArray(od->xbuffer, 3, &od->xLastMax, &ix);
                // od->tLastMax = od->tbuffer[ix];
            }

            // //local min check in fVarIx - one between each max - simply overwrite tLastMin, xLastMin
            // if (od->dxbuffer[1] <= op->eps_dx && od->dxbuffer[2] > op->eps_dx)
            // {
            // }

            //Check for downward threshold crossing (end of "active phase")
            if (xi[op->fVarIx] <= od->xDown && dxi[op->fVarIx] >= od->dxDown && od->inUpstate)
            {
                //simplest: time and state of trajectory point that landed in the neighborhood.
                // od->tThisDown=*ti;

                //alts: linearly interp to get more accurate tThisUp within the last time step: xi-1, xi, ti-1, ti; get t @ x=xUp
                od->tThisDown = linearInterp(od->xbuffer[2], xi[op->fVarIx], od->tbuffer[2], *ti, od->xDown);
                od->inUpstate = false;
            }
        }

    }

    //advance solution buffer last (so really have: buffer[0,1,2],xi) TODO: alt - update buffers immediately after step, so only use buffers for event/obs? Related: FSAL.
    od->tbuffer[0] = od->tbuffer[1];
    od->tbuffer[1] = od->tbuffer[2];
    od->tbuffer[2] = *ti;
    od->xbuffer[0] = od->xbuffer[1];
    od->xbuffer[1] = od->xbuffer[2];
    od->xbuffer[2] = xi[op->fVarIx];
    od->dxbuffer[0] = od->dxbuffer[1];
    od->dxbuffer[1] = od->dxbuffer[2];
    od->dxbuffer[2] = dxi[op->fVarIx];
    od->maxDt = fmax(od->tbuffer[2]-od->tbuffer[1],od->maxDt);
    od->minDt = fmin(od->tbuffer[2]-od->tbuffer[1],od->minDt);

    if (eventOccurred)
    { //reset event-based quantities
    }
}

//Perform and post-integration cleanup and write desired features into the global array F
void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype *F, int i, int nPts)
{
    //Number of features is determined by this function. Must hardcode that number into the host program in order to allocate memory for F...
    int ix = 0;
    F[ix++ * nPts + i] = od->xGlobalMax - od->xGlobalMin;
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->period[0] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->period[1] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->period[2] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->nMaxima[0] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->nMaxima[1] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->nMaxima[2] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->upDuration[0] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->upDuration[1] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->upDuration[2] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->downDuration[0] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->downDuration[1] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->downDuration[2] : RCONST(0.0);
    F[ix++ * nPts + i] = od->dxGlobalMax;
    F[ix++ * nPts + i] = od->dxGlobalMin;
    F[ix++ * nPts + i] = od->xTrajectoryMean;
    F[ix++ * nPts + i] = od->eventcount;
    F[ix++ * nPts + i] = od->stepcount;
    F[ix++ * nPts + i] = od->maxDt;
    F[ix++ * nPts + i] = od->minDt;
}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype *tspan)
{
    // od->stepcount = 0;
    //shift all time-based observer members left by [tf-t0]
    realtype T = *ti - tspan[0];
    od->tLastEvent -= T; //i.e. tLastUp
    od->tThisDown -= T;
    // od->tLastMax-=T;
    // od->tLastMin-=T;
    for (int j = 0; j < 3; ++j)
    {
        od->tbuffer[j] -= T;
    }
    // od->maxDt = -BIG_REAL;
    // od->minDt = BIG_REAL;
}

#endif // OBSERVER_THRESHOLD_2_H
