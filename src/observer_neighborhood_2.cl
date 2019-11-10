//

#ifndef OBSERVER_NEIGHBORHOOD_2_H
#define OBSERVER_NEIGHBORHOOD_2_H

#include "clODE_utilities.cl"
#include "observers.cl"
#include "realtype.cl"

#define TWO_PASS_EVENT_DETECTOR

//typedef to select the correct observer struct from the header
typedef struct ObserverData
{

    realtype tbuffer[3];
    realtype xbuffer[3 * N_VAR];
    realtype dxbuffer[3];

    realtype x0[N_VAR]; //point for neighborhood return
    // realtype xLast[N_VAR];
    // realtype xMaxDelta[N_VAR]; //store the maximum delta-x = xi[j]-xLast[j] per time step in each dimension (for size of epsilon "ball")

    realtype xTrajectoryMax[N_VAR];
    realtype xTrajectoryMin[N_VAR];
    realtype xTrajectoryMean[N_VAR];
    realtype dxTrajectoryMax[N_VAR];
    realtype dxTrajectoryMin[N_VAR];
    realtype xRange[N_VAR]; //to store warmup
    // realtype dxRange[N_VAR];

    realtype auxTrajectoryMax[N_AUX];
    realtype auxTrajectoryMin[N_AUX];
    realtype auxTrajectoryMean[N_AUX];

    //period-wise features
    realtype nMaxima[3]; //max/min/mean
    realtype period[3];  //max/min/mean
    // realtype IMI[3]; //max/min/mean
    // realtype amp[3]; //max/min/mean
    // realtype xMax[3]; //max/min/mean
    // realtype xMin[3]; //max/min/mean
    realtype stepDt[3]; //max/min/mean

    realtype tLastEvent;
    realtype xThreshold;
    realtype thisNormXdiff;
    realtype lastNormXdiff;

    // realtype xLastMax;
    // realtype tLastMax;

    // realtype xLastMin;
    // realtype tLastMin;

    int thisNMaxima;
    int eventcount;
    int stepcount;
    bool buffer_filled;
    bool foundX0;
    bool isInNhood;

} ObserverData;
// ObserverData size: 3*int + (4*3 + 7*nVar + 9)*realtype + 3*bool

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    od->eventcount = 0;
    od->stepcount = 0;
    od->foundX0 = false;
    od->isInNhood = false;

    //put x0 in the leading position of the solution buffer
    od->buffer_filled = false;
    od->tbuffer[2] = *ti;
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xbuffer[j * 3 + 2] = xi[j];
    }
    od->dxbuffer[2] = dxi[op->fVarIx];

    for (int j = 0; j < N_VAR; ++j)
    {
        // od->x0[j] = RCONST(0.0); //not needed - set first time anyway.
        // od->x0[j] = xi[j];
        od->xTrajectoryMean[j] = RCONST(0.0); //not needed - set first time anyway.
        od->xTrajectoryMax[j] = -BIG_REAL;
        od->xTrajectoryMin[j] = BIG_REAL;
        od->dxTrajectoryMax[j] = -BIG_REAL;
        od->dxTrajectoryMin[j] = BIG_REAL;
    }
    for (int j = 0; j < N_AUX; ++j)
    {
        od->auxTrajectoryMax[j] = -BIG_REAL;
        od->auxTrajectoryMin[j] = BIG_REAL;
        od->auxTrajectoryMean[j] = RCONST(0.0); //not needed - set first time anyway.
    }
    // od->tLastEvent=RCONST(0.0);

    // od->xLastMax = -BIG_REAL;
    // od->tLastMax = RCONST(0.0);

    // od->xLastMin = BIG_REAL;
    // od->tLastMin = RCONST(0.0);

    od->thisNMaxima = 0;

    od->nMaxima[0] = -BIG_REAL;
    od->nMaxima[1] = BIG_REAL;
    od->nMaxima[2] = RCONST(0.0);

    od->period[0] = -BIG_REAL;
    od->period[1] = BIG_REAL;
    od->period[2] = RCONST(0.0);

    // od->IMI[0]=-BIG_REAL;
    // od->IMI[1]= BIG_REAL;
    // od->IMI[2]= RCONST(0.0);

    // od->amp[0]=-BIG_REAL;
    // od->amp[1]= BIG_REAL;
    // od->amp[2]= RCONST(0.0);

    // od->xMax[0]=-BIG_REAL;
    // od->xMax[1]= BIG_REAL;
    // od->xMax[2]= RCONST(0.0);

    // od->xMin[0]=-BIG_REAL;
    // od->xMin[1]= BIG_REAL;
    // od->xMin[2]= RCONST(0.0);

    od->stepDt[0] = -BIG_REAL;
    od->stepDt[1] = BIG_REAL;
    od->stepDt[2] = RCONST(0.0);

    // od->xTrajectoryMean = RCONST(0.0);
    // od->xGlobalMax = -BIG_REAL;
    // od->xGlobalMin = BIG_REAL;
    // od->dxGlobalMax = -BIG_REAL;
    // od->dxGlobalMin = BIG_REAL;
}

//restricted per-timestep update of observer data for initializing event detector
// - get extent of trajectory in state space, and max/min slopes
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    for (int j = 0; j < N_VAR; ++j)
    {
        // od->xMaxDelta[j]=fmax( fabs( xi[j] - od->xLast[j] ), od->xMaxDelta[j] );
        // od->xLast[j]=xi[j];

        od->xTrajectoryMax[j] = fmax(xi[j], od->xTrajectoryMax[j]);
        od->xTrajectoryMin[j] = fmin(xi[j], od->xTrajectoryMin[j]);
        // od->dxTrajectoryMax[j]=fmax(dxi[j], od->dxTrajectoryMax[j]);
        // od->dxTrajectoryMin[j]=fmin(dxi[j], od->dxTrajectoryMin[j]);
    }

    // //find abs min as x0
    // if (xi[op->eVarIx] < od->x0[op->eVarIx]) {
    // 	for (int j=0; j<N_VAR;++j) {
    // 		od->x0[j]=xi[j];
    // 	}
    // 	od->tLastEvent=*ti;
    // }
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{

    for (int j = 0; j < N_VAR; ++j)
    {
        // od->xMaxDelta[j]=RCONST(0.0); //maximum magnitude of step taken by solver in each state variable direction
        // od->x0[j]=xi[j];
        // od->xLast[j]=xi[j];
        od->xRange[j] = od->xTrajectoryMax[j] - od->xTrajectoryMin[j];
        // od->dxRange[j]=od->dxTrajectoryMax[j]-od->dxTrajectoryMin[j];
    }

    od->xThreshold = od->xTrajectoryMin[op->eVarIx] + op->xDownThresh * od->xRange[op->eVarIx]; 
    //TODO: add downward threshold too, for "up" and "down" state durations
    // od->xThreshold=od->xTrajectoryMin[op->eVarIx] + op->xUpThresh*od->xRange[op->eVarIx];

    // od->isInNhood = true; //we are at the center of the neighborhood, x0.
    // od->tLastEvent=*ti;

    for (int j = 0; j < N_VAR; ++j)
    {
        od->xTrajectoryMax[j] = -BIG_REAL;
        od->xTrajectoryMin[j] = BIG_REAL;
        // od->dxTrajectoryMax[j] = -BIG_REAL;
        // od->dxTrajectoryMin[j] = BIG_REAL;
    }
    // for (int j = 0; j < N_AUX; ++j)
    // {
    //     od->auxTrajectoryMax[j] = -BIG_REAL;
    //     od->auxTrajectoryMin[j] = BIG_REAL;
    // }
}

//check for entry into "epsilon ball" surrounding od->x0
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    //minimum amplitude check in variable: fVarIx 
    if (od->buffer_filled == 1 && od->xRange[op->fVarIx] > op->minXamp)
    {
        if (od->foundX0)
        {
            //save values for comparison from last timepoint
            bool lastInNhood = od->isInNhood;
            od->lastNormXdiff = od->thisNormXdiff;

            realtype thisXdiff[N_VAR];
            for (int j = 0; j < N_VAR; ++j)
                thisXdiff[j] = fabs(xi[j] - od->x0[j]) / od->xRange[j];

            // od->thisNormXdiff=norm_1(thisXdiff, N_VAR);
            od->thisNormXdiff = norm_2(thisXdiff, N_VAR);
            od->isInNhood = od->thisNormXdiff <= op->nHoodRadius; //L-2 norm ball
            return (od->isInNhood & ~lastInNhood); //event only on entry
        }
    }
    return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
// - get (t,x) at local max, compute: IMI, tMaxMin, amp
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    realtype tThisEvent;

    ++od->eventcount;

    if (od->eventcount == 1)
    {
        tThisEvent = *ti;
    }
    else if (od->eventcount > 1)
    {

        //simplest: time and state of trajectory point that landed in the neighborhood.
        // tThisEvent=*ti;

        //alts: linearly interp to get tThisEvent using thisNormXDiff/lastNormXDiff; use norm values as x, since we know the interp value of that
        tThisEvent = linearInterp(od->lastNormXdiff, od->thisNormXdiff, od->tbuffer[2], *ti, op->nHoodRadius);

        od->nMaxima[0] = fmax((realtype)od->thisNMaxima, od->nMaxima[0]); //cast to realtype (for mean)
        od->nMaxima[1] = fmin((realtype)od->thisNMaxima, od->nMaxima[1]);
        runningMean(&od->nMaxima[2], (realtype)od->thisNMaxima, od->eventcount - 1);

        realtype thisPeriod = tThisEvent - od->tLastEvent;
        //realtype thisPeriod=*ti-tThisEvent;
        od->period[0] = fmax(thisPeriod, od->period[0]);
        od->period[1] = fmin(thisPeriod, od->period[1]);
        runningMean(&od->period[2], thisPeriod, od->eventcount - 1);
    }

    od->tLastEvent = tThisEvent;
    od->thisNMaxima = 0;

    return false; //not terminal
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
// - advance solution/slope buffers
// - check for any intermediate special points & store their info
// - reset intermediates upon local max event detection
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred)
{

    // realtype tThisMax, xThisMax;
    // realtype thisIMI, thisTMaxMin, thisAmp;

    // ++od->stepcount;

    //global extent of all vars, var slopes, and aux vars
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xTrajectoryMax[j] = fmax(xi[j], od->xTrajectoryMax[j]);
        od->xTrajectoryMin[j] = fmin(xi[j], od->xTrajectoryMin[j]);
        runningMean(&od->xTrajectoryMean[j], xi[j], od->stepcount);
        od->dxTrajectoryMax[j] = fmax(dxi[j], od->dxTrajectoryMax[j]);
        od->dxTrajectoryMin[j] = fmin(dxi[j], od->dxTrajectoryMin[j]);
    }
    for (int j = 0; j < N_AUX; ++j)
    {
        od->auxTrajectoryMax[j] = fmax(auxi[j], od->auxTrajectoryMax[j]);
        od->auxTrajectoryMin[j] = fmin(auxi[j], od->auxTrajectoryMin[j]);
        runningMean(&od->auxTrajectoryMean[j], auxi[j], od->stepcount);
    }

    if (od->stepcount > 2 && od->buffer_filled == 0)
    {
        od->buffer_filled = 1;
    }

    if (od->buffer_filled == 1)
    {
        //check for x0
        if (!od->foundX0)
        {   //x0 is first time dropping below threshold in x[eVarIx] - xbuffer holds previous xi
            float lastX=od->xbuffer[op->eVarIx * 3 + 2], thisX=xi[op->eVarIx];
            if (lastX > od->xThreshold && thisX < od->xThreshold)
            {
                od->foundX0 = true;
                od->isInNhood = true;

                //record x0
                for (int j = 0; j < N_VAR; ++j)
                    od->x0[j] = xi[j];

                //get plane equation with point X0=Xi and normal vector (X0 - Xi-1)
            }
        }

        //local max check in fVarIx - one between each min - simply overwrite tLastMax, xLastMax
        // if (od->dxbuffer[1] > -op->eps_dx && od->dxbuffer[2] < -op->eps_dx)
        if (od->dxbuffer[2] > -op->eps_dx && dxi[op->fVarIx] < -op->eps_dx)
        {
            od->thisNMaxima++;
            // int ix;
            // realtype thisXbuffer[N_VAR];
            // for (int j = 0; j < 3; ++j)
            // 	thisXbuffer[j] = od->xbuffer[op->fVarIx * 3 + j];
            // maxOfArray(thisXbuffer, 3, &od->xLastMax, &ix);
            // od->tLastMax = od->tbuffer[ix];

            // if (od->thisNMaxima>1)
            // {

            // }
        }

        // //local min check - one between each max - simply overwrite tLastMin, xLastMin
        // if (od->dxbuffer[1] <= op->eps_dx && od->dxbuffer[2] > op->eps_dx ) {
        // 	int ix;
        // 	minOfArray(od->xbuffer, 3, &od->xLastMin, &ix);
        // 	od->tLastMin=od->tbuffer[ix];
        // }

        // }
        // else
        // {

        // 	for (int j=0; j<N_VAR; ++j){
        // 	// 	od->xMaxDelta[j]=od->xLast[j]-xi[j];//fmax( fabs(od->xLast[j]-xi[j]), od->xMaxDelta[j] );
        // 		// od->x0[j]=xi[j];
        // 		od->xLast[j]=xi[j];
        // 		// od->tLastEvent=*ti;
        // 	}
    }

    //advance solution buffer
    od->tbuffer[0] = od->tbuffer[1];
    od->tbuffer[1] = od->tbuffer[2];
    od->tbuffer[2] = *ti;
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xbuffer[j * 3 + 0] = od->xbuffer[j * 3 + 1];
        od->xbuffer[j * 3 + 1] = od->xbuffer[j * 3 + 2];
        od->xbuffer[j * 3 + 2] = xi[j];
    }
    od->dxbuffer[0] = od->dxbuffer[1];
    od->dxbuffer[1] = od->dxbuffer[2];
    od->dxbuffer[2] = dxi[op->fVarIx];

    //record actual dt
    realtype thisDt = od->tbuffer[2] - od->tbuffer[1];
    od->stepDt[0] = fmax(thisDt, od->stepDt[0]);
    od->stepDt[1] = fmin(thisDt, od->stepDt[1]);
    runningMean(&od->stepDt[2], thisDt, od->stepcount);
    
}

//Perform and post-integration cleanup and write desired features into the global array F
void finalizeFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __global realtype *F, int i, int nPts)
{
    int ix = 0;
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->period[0] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->period[1] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->period[2] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->nMaxima[0] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->nMaxima[1] : RCONST(0.0);
    F[ix++ * nPts + i] = od->eventcount > 1 ? od->nMaxima[2] : RCONST(0.0);
    for (int j = 0; j < N_VAR; ++j) //5*N_VAR
    {
        F[ix++ * nPts + i] = od->xTrajectoryMax[j];
        F[ix++ * nPts + i] = od->xTrajectoryMin[j];
        F[ix++ * nPts + i] = od->xTrajectoryMean[j];
        F[ix++ * nPts + i] = od->dxTrajectoryMax[j];
        F[ix++ * nPts + i] = od->dxTrajectoryMin[j];
    }
    for (int j = 0; j < N_AUX; ++j) //3*N_AUX
    {
        F[ix++ * nPts + i] = od->auxTrajectoryMax[j];
        F[ix++ * nPts + i] = od->auxTrajectoryMin[j];
        F[ix++ * nPts + i] = od->auxTrajectoryMean[j];
    }
    F[ix++ * nPts + i] = od->eventcount;
    F[ix++ * nPts + i] = od->stepcount;
    F[ix++ * nPts + i] = od->stepDt[0];
    F[ix++ * nPts + i] = od->stepDt[1];
    F[ix++ * nPts + i] = od->stepDt[2];
}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype *tspan)
{
    //shift all time-based observer members left by [tf-t0]
    realtype T = *ti - tspan[0];
    od->tLastEvent -= T;
    // od->tLastMax-=T;
    // od->tLastMin-=T;
    for (int j = 0; j < 3; ++j)
    {
        od->tbuffer[j] -= T;
    }
}

#endif // OBSERVER_NEIGHBORHOOD_2_H
