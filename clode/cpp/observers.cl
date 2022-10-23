#ifndef OBSERVERS_H_
#define OBSERVERS_H_

/* "Observer" measures features of the ODE solution as it is being integrated
 * The observer consists of a data structure and several functions: 
 * 
 * - initializeObserverData: set up the data structure to sensible values
 * - warmupObserverData: for two-pass event detectors - restricted data collection about trajectory during a first pass ODE solve
 * - updateObserverData: per-timestep update of data structure
 * - initializeEventDetector: set any values needed to do selected type of event detection (possibly using warmup data)
 * - eventFunction: check for an event. Optionally refine location of event within timestep. Compute event-based quantities
 * - computeEventFeatures: when event is detected, compute desired per-event features
 * - finalizeFeatures: post-integration cleanup and write to global feature array
 */

#include "realtype.cl"

#ifdef __cplusplus
template <typename realtype>
#endif
struct ObserverParams
{
    int eVarIx; //variable for event detection
    int fVarIx; //variable for features

    int maxEventCount; //time loop limiter
    realtype minXamp;  //consider oscillations lower than this to be steady state (return mean X)
    realtype minIMI;

    //neighborhood return map
    realtype nHoodRadius;

    //section. Two interpretations: absolute, relative
    realtype xUpThresh;
    realtype xDownThresh;
    realtype dxUpThresh;
    realtype dxDownThresh;

    //local extremum - tolerance for zero crossing of dx - for single precision: if RHS involves sum of terms of O(1), dx=zero is noise at O(1e-7)
    realtype eps_dx;
};


#ifdef __cplusplus
//info about available observers for access in C++
typedef struct ObserverInfo 
{
	std::string define;
	size_t observerDataSizeFloat;
	size_t observerDataSizeDouble;
	std::vector<std::string> featureNames;
} ObserverInfo;

#endif


////////////////////////////////////////////////
// one-pass detectors
////////////////////////////////////////////////

//basic detectors: no events. Measure extent of state space explored, max/min/mean x and aux, max/min dx
#include "observers/observer_basic.clh" //one variable, specified by fVarIx
#include "observers/observer_basic_allVar.clh"
#include "observers/observer_local_maximum.clh"

// #include "observers/observer_threshold_1.clh" //not implemented

//Event is the return of trajectory to small neighborhood of Xstart (specified state, eg. x0)
// to select Xstart: local min (e.g. of user-selected slow variable) 
#include "observers/observer_neighborhood_1.clh"

////////////////////////////////////////////////
// two-pass detectors
////////////////////////////////////////////////

//Threshold-based event detection with relative thresholds in a specified variable xi.
// First pass to measure the extent of state-space trajectory visits, then compute thresholds as fractions of range
#include "observers/observer_threshold_2.clh"

// // First pass to measure the extent of state-space trajectory visits, decide x0 based on threshold. check each crossing, period when norm(x-x0)<tol
// #include "observers/observer_poincare_2.clh"

//Use a first pass to find a good Xstart (e.g. absolute drop below 0.5*range of slowest variable)
#include "observers/observer_neighborhood_2.clh"



// collect available methods into "name"-ObserverInfo map, for C++ side access. Must come after including all the getObserverInfo_functions.
#ifdef __cplusplus
static void getObserverDefineMap(const ProblemInfo pi, const int fVarIx, const int eVarIx, std::map<std::string, ObserverInfo> &observerDefineMap, std::vector<std::string> &availableObserverNames) 
{

    std::map<std::string, ObserverInfo> newMap;
    newMap["basic"]=getObserverInfo_basic(pi, fVarIx, eVarIx);
    newMap["basicall"]=getObserverInfo_basicAll(pi, fVarIx, eVarIx);
    newMap["localmax"]=getObserverInfo_localmax(pi, fVarIx, eVarIx);
    newMap["nhood1"]=getObserverInfo_nhood1(pi, fVarIx, eVarIx);
    newMap["nhood2"]=getObserverInfo_nhood2(pi, fVarIx, eVarIx);
    newMap["thresh2"]=getObserverInfo_thresh2(pi, fVarIx, eVarIx);

//export vector of names for access in C++
std::vector<std::string> newNames;
for (auto const& element : newMap)
    newNames.push_back(element.first);

observerDefineMap=newMap;
availableObserverNames=newNames;
}
#endif


#endif //OBSERVERS_H_

/*
struct SolBuffer {
	realtype t[BUFFER_SIZE];
	realtype x[BUFFER_SIZE][N_VAR];
	realtype dx[BUFFER_SIZE][N_VAR];
	realtype aux[BUFFER_SIZE][N_AUX];
};

void updateSolutionBuffer(struct SolBuffer *sb, realtype *ti, realtype xi[], realtype dxi[], realtype auxi[]) {
	for (int i=0; i<BUFFER_SIZE-1; ++i) {
		sb->t[i]=sb->t[i+1];
		for (int j=0; j<N_VAR; ++j) {
			sb->x[i][j]=xi[j];
			sb->dx[i][j]=dxi[j];
		}
		for (int j=0; j<N_AUX; ++j) {
			sb->aux[i][j]=auxi[j];
		}
	}
	sb->t[BUFFER_SIZE]=*ti;
	for (int j=0; j<N_VAR; ++j) {
		sb->x[BUFFER_SIZE][j]=xi[j];
		sb->dx[BUFFER_SIZE][j]=dxi[j];
	}
	for (int j=0; j<N_AUX; ++j) {
		sb->aux[BUFFER_SIZE][j]=auxi[j];
	}
}
*/
