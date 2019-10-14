//

#ifndef OBSERVER_NEIGHBORHOOD_2_H
#define OBSERVER_NEIGHBORHOOD_2_H

#include "clODE_utilities.cl"
#include "observers.cl"

// #define TWO_PASS_EVENT_DETECTOR

//broken attempt to use x0.

//typedef to select the correct observer struct from the header
typedef struct ObserverData
{

    realtype tbuffer[3];
    realtype xbuffer[3 * N_VAR];
    realtype dxbuffer[3];

    realtype x0[N_VAR]; //point for neighborhood return
    realtype xLast[N_VAR];
    realtype xMaxDelta[N_VAR]; //store the maximum delta-x = xi[j]-xLast[j] per time step in each dimension (for size of epsilon "ball")
    realtype tLastX0;

    // realtype xTrajectoryMax[N_VAR];
    // realtype xTrajectoryMin[N_VAR];
    // realtype xRange[N_VAR];
    // realtype dxTrajectoryMax[N_VAR];
    // realtype dxTrajectoryMin[N_VAR];
    // realtype dxRange[N_VAR];

    realtype nMaxima[3]; //max/min/mean
    realtype period[3];  //max/min/mean
    // realtype IMI[3]; //max/min/mean
    // realtype amp[3]; //max/min/mean
    // realtype xMax[3]; //max/min/mean
    // realtype xMin[3]; //max/min/mean

    realtype dxMax;
    realtype dxMin;
    realtype xTrajectoryMean;

    realtype xLastMax;
    realtype tLastMax;

    realtype xLastMin;
    realtype tLastMin;

    int thisNMaxima;
    int eventcount;
    int stepcount;
    bool buffer_filled;
    bool isInNhood;

} ObserverData;
// ObserverData size: 3*int + (3*3 + 3*nVar + 8 + 2*3)*realtype + 2*bool

//set initial values to relevant fields in ObserverData
void initializeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{

    od->eventcount = 0;
    od->stepcount = 0;

    //put x0 in the leading position of the solution buffer
    od->buffer_filled = false;
    od->tbuffer[2] = *ti;
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xbuffer[j * 3 + 2] = xi[j];
    }
    od->dxbuffer[2] = dxi[op->fVarIx];

    // for (int j=0; j<N_VAR;++j) {
    // 	od->x0[j]=xi[j];
    // 	od->xLast[j]=xi[j];
    // 	od->xMaxDelta[j]=RCONST(0.0); //maximum magnitude of step taken by solver in each state variable direction

    // 	// od->xTrajectoryMax[j]=-BIG_REAL;
    // 	// od->xTrajectoryMin[j]= BIG_REAL;
    // 	// od->dxTrajectoryMax[j]=-BIG_REAL;
    // 	// od->dxTrajectoryMin[j]= BIG_REAL;
    // }
    // od->tLastX0=RCONST(0.0);

    od->xLastMax = -BIG_REAL;
    od->tLastMax = RCONST(0.0);

    od->xLastMin = BIG_REAL;
    od->tLastMin = RCONST(0.0);

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

    od->xTrajectoryMean = RCONST(0.0);
    od->dxMax = -BIG_REAL;
    od->dxMin = BIG_REAL;
}

//restricted per-timestep update of observer data for initializing event detector
// - get extent of trajectory in state space, and max/min slopes
void warmupObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{
    for (int j = 0; j < N_VAR; ++j)
    {
        od->xMaxDelta[j] = MAX(fabs(xi[j] - od->xLast[j]), od->xMaxDelta[j]);
        od->xLast[j] = xi[j];

        // od->xTrajectoryMax[j]=MAX(xi[j], od->xTrajectoryMax[j]);
        // od->xTrajectoryMin[j]=MIN(xi[j], od->xTrajectoryMin[j]);
        // od->dxTrajectoryMax[j]=MAX(dxi[j], od->dxTrajectoryMax[j]);
        // od->dxTrajectoryMin[j]=MIN(dxi[j], od->dxTrajectoryMin[j]);
    }

    //find abs min as x0
    // if (xi[op->fVarIx]<od->x0[op->fVarIx]) {
    // 	for (int j=0; j<N_VAR;++j) {
    // 		od->x0[j]=xi[j];
    // 	}
    // }
}

//process warmup data to compute relevant event detector quantities (e.g. thresholds)
void initializeEventDetector(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{

    for (int j = 0; j < N_VAR; ++j)
    {
        od->xMaxDelta[j] = RCONST(0.0); //maximum magnitude of step taken by solver in each state variable direction
        od->x0[j] = xi[j];
        od->xLast[j] = xi[j];
        // od->xRange[j]=od->xTrajectoryMax[j]-od->xTrajectoryMin[j];
        // od->dxRange[j]=od->dxTrajectoryMax[j]-od->dxTrajectoryMin[j];
    }
    od->isInNhood = true; //we are at the center of the neighborhood, x0.
    od->tLastX0 = *ti;
}

//check for entry into "epsilon ball" surrounding od->x0
bool eventFunction(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{

    if (od->buffer_filled == 1)
    {
        bool lastWasInNhood = od->isInNhood;

        realtype thisXdiff[N_VAR];
        // od->isInNhood=true;
        for (int j = 0; j < N_VAR; ++j)
        {
            // od->isInNhood=od->isInNhood & ( fabs(od->x0[j]-xi[j]) < RCONST(10.0)*od->xMaxDelta[j] );
            thisXdiff[j] = fabs(od->x0[j] - xi[j]) / od->xMaxDelta[j];
        }

        // od->isInNhood=norm_1( thisXdiff , N_VAR ) < RCONST(5.0); //current point is within an L-1 norm ball of x0
        od->isInNhood = norm_2(thisXdiff, N_VAR) < RCONST(4.0); //L-2 norm ball

        return (od->isInNhood & ~lastWasInNhood); //event only on entry
    }
    else //collect the larger of x deltas near x0
    {
        for (int j = 0; j < N_VAR; ++j)
        {
            od->xMaxDelta[j] = MAX(fabs(od->xLast[j] - xi[j]), od->xMaxDelta[j]);
            // od->xLast[j]=xi[j];
        }
    }

    // if (od->stepcount>2){return (true);}

    return false;
}

//When an event is detected, computes desired event-based features. returns true if a terminal event was reached
// - get (t,x) at local max, compute: IMI, tMaxMin, amp
bool computeEventFeatures(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op)
{

    ++od->eventcount;

    if (od->eventcount > 1)
    {

        od->nMaxima[0] = MAX((realtype)od->thisNMaxima, od->nMaxima[0]); //cast to realtype (for mean)
        od->nMaxima[1] = MIN((realtype)od->thisNMaxima, od->nMaxima[1]);
        runningMean(&od->nMaxima[2], (realtype)od->thisNMaxima, od->eventcount - 1);

        float thisPeriod = *ti - od->tLastX0;
        od->period[0] = MAX(thisPeriod, od->period[0]);
        od->period[1] = MIN(thisPeriod, od->period[1]);
        runningMean(&od->period[2], thisPeriod, od->eventcount - 1);
    }

    return false; //not terminal
}

//full per-timestep update of observer data. If an event occurred this timestep, event-based observer data is reset. Per-timestep features are computed here.
// - advance solution/slope buffers
// - check for any intermediate special points & store their info
// - reset intermediates upon local max event detection
void updateObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, bool eventOccurred)
{

    //reset intermediate feature storage
    if (eventOccurred)
    {
        od->tLastX0 = *ti;
        od->thisNMaxima = 0;
    }

    // realtype tThisMax, xThisMax;
    // realtype thisIMI, thisTMaxMin, thisAmp;

    // ++od->stepcount;
    if (od->stepcount > 2 && od->buffer_filled == 0)
    {
        od->buffer_filled = 1;
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

    //global dxMax, dxMin
    od->dxMax = MAX(od->dxMax, dxi[op->fVarIx]);
    od->dxMin = MIN(od->dxMin, dxi[op->fVarIx]);
    runningMean(&od->xTrajectoryMean, xi[op->fVarIx], od->stepcount);

    if (od->buffer_filled == 1)
    {

        //local max check - one between each min - simply overwrite tLastMax, xLastMax
        if (od->dxbuffer[1] >= -op->eps_dx && od->dxbuffer[2] < -op->eps_dx)
        {
            od->thisNMaxima++;
            int ix;
            float thisXbuffer[N_VAR];
            for (int j = 0; j < 3; ++j)
                thisXbuffer[j] = od->xbuffer[op->fVarIx * 3 + j];
            maxOfArray(thisXbuffer, 3, &od->xLastMax, &ix);
            od->tLastMax = od->tbuffer[ix];

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
    }
    else
    {

        for (int j = 0; j < N_VAR; ++j)
        {
            // 	od->xMaxDelta[j]=od->xLast[j]-xi[j];//MAX( fabs(od->xLast[j]-xi[j]), od->xMaxDelta[j] );
            // od->x0[j]=xi[j];
            od->xLast[j] = xi[j];
            // od->tLastX0=*ti;
        }
    }
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
    // F[ix++*nPts+i]=od->eventcount>0?od->xMax[0]:xi[op->fVarIx];
    // F[ix++*nPts+i]=od->eventcount>0?od->xMax[1]:xi[op->fVarIx];
    // F[ix++*nPts+i]=od->eventcount>0?od->xMax[2]:xi[op->fVarIx];
    // F[ix++*nPts+i]=od->eventcount>0?od->xMin[0]:xi[op->fVarIx];
    // F[ix++*nPts+i]=od->eventcount>0?od->xMin[1]:xi[op->fVarIx];
    // F[ix++*nPts+i]=od->eventcount>0?od->xMin[2]:xi[op->fVarIx];

    // F[ix++*nPts+i]=N_VAR;
    // int n=2;
    // F[ix++*nPts+i]=od->xbuffer[0*3+n];
    // F[ix++*nPts+i]=od->xbuffer[1*3+n];
    // F[ix++*nPts+i]=od->xbuffer[2*3+n];
    // F[ix++*nPts+i]=od->xbuffer[3*3+n];
    // F[ix++*nPts+i]=od->xbuffer[0];
    // F[ix++*nPts+i]=od->xbuffer[3];
    // F[ix++*nPts+i]=od->xbuffer[6];
    // F[ix++*nPts+i]=od->xbuffer[9];
    // F[ix++*nPts+i]=xi[0];
    // F[ix++*nPts+i]=xi[1];
    // F[ix++*nPts+i]=xi[2];
    // F[ix++*nPts+i]=xi[3];
    // F[ix++*nPts+i]=od->x0[0];
    // F[ix++*nPts+i]=od->x0[1];
    // F[ix++*nPts+i]=od->x0[2];
    // F[ix++*nPts+i]=od->x0[3];
    // F[ix++*nPts+i]=od->xLast[0];
    // F[ix++*nPts+i]=od->xLast[1];
    // F[ix++*nPts+i]=od->xLast[2];
    // F[ix++*nPts+i]=od->xLast[3];
    // F[ix++*nPts+i]=od->xMaxDelta[0];
    // F[ix++*nPts+i]=od->xMaxDelta[1];
    // F[ix++*nPts+i]=od->xMaxDelta[2];
    // F[ix++*nPts+i]=od->xMaxDelta[3];
    // F[ix++*nPts+i]=fabs(od->x0[0]-xi[0]);
    // F[ix++*nPts+i]=fabs(od->x0[1]-xi[1]);
    // F[ix++*nPts+i]=fabs(od->x0[2]-xi[2]);
    // F[ix++*nPts+i]=fabs(od->x0[3]-xi[3]);
    // F[ix++*nPts+i]=fabs(od->xLast[0]-xi[0]);
    // F[ix++*nPts+i]=fabs(od->xLast[1]-xi[1]);
    // F[ix++*nPts+i]=fabs(od->xLast[2]-xi[2]);
    // F[ix++*nPts+i]=fabs(od->xLast[3]-xi[3]);
    // F[ix++*nPts+i]=fabs(od->x0[0]-xi[0])/od->xMaxDelta[0];
    // F[ix++*nPts+i]=fabs(od->x0[1]-xi[1])/od->xMaxDelta[1];
    // F[ix++*nPts+i]=fabs(od->x0[2]-xi[2])/od->xMaxDelta[2];
    // F[ix++*nPts+i]=fabs(od->x0[3]-xi[3])/od->xMaxDelta[3];
    F[ix++ * nPts + i] = od->dxMax;
    F[ix++ * nPts + i] = od->dxMin;
    F[ix++ * nPts + i] = od->xTrajectoryMean;
    F[ix++ * nPts + i] = od->eventcount;
    F[ix++ * nPts + i] = od->stepcount;
}

//Perform and post-integration cleanup of observer data to ensure it is ready for continuation if needed
void finalizeObserverData(realtype *ti, realtype xi[], realtype dxi[], realtype auxi[], ObserverData *od, __constant struct ObserverParams *op, __constant realtype tspan[])
{
    //shift all time-based observer members left by [tf-t0]
    realtype T = *ti - tspan[0];
    od->tLastMax = od->tLastMax - T;
    od->tLastMin = od->tLastMin - T;
    for (int j = 0; j < 3; ++j)
    {
        od->tbuffer[j] = od->tbuffer[j] - T;
    }
}

#endif // OBSERVER_NEIGHBORHOOD_2_H
