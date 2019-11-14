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
	size_t observerDataSize;
	std::vector<std::string> featureNames;
} ObserverInfo;

#endif


////////////////////////////////////////////////
// one-pass detectors
////////////////////////////////////////////////

//basic detector: no events. All features are computed per step (max/min/mean vars/aux)
#include "observers/observer_basic.clh"

#ifdef OBSERVER_BASIC_ALLVAR
#include "observers/observer_basic_allVar.clh"
#endif

//Event is the detection of local extremum (max and/or min) in a specified variable xi
#ifdef OBSERVER_LOCAL_MAX
#include "observers/observer_local_maximum.clh"
#endif

//Threshold-based event detection with absolute thresholds in a specified variable xi. Three flavors:
// 1) xup = value of xi. (simple threshold)
// 2) xup = value of xi, xdown = value of xi < xup (Shmitt trigger)
// 3) xup, xdown, dxup, dxdown. (Shmitt with slope thresholds)
#ifdef OBSERVER_THRESHOLD_1
#include "observers/observer_threshold_1.clh"
#endif

//Event is the return of trajectory to small neighborhood of Xstart (specified state, eg. x0)
// can use transient to approach periodic solution, then use x0
#ifdef OBSERVER_NEIGHBORHOOD_1
#include "observers/observer_neighborhood_1.clh"
#endif

////////////////////////////////////////////////
// two-pass detectors
////////////////////////////////////////////////

//TODO: make a separate kernel to use for warmup (to break the computation into parts so host doesn't freeze)

//Threshold-based event detection with relative thresholds in a specified variable xi.
// Need to measure the extent of state-space trajectory visits, then compute thresholds as fractions of range
#ifdef OBSERVER_THRESHOLD_2
#include "observers/observer_threshold_2.clh"
#endif

//Use a first pass to find a good Xstart (e.g. absolute min of slowest variable)
#ifdef OBSERVER_NEIGHBORHOOD_2
#include "observers/observer_neighborhood_2.clh"
#endif



// collect available methods into "name"-ObserverInfo map, for C++ side access. Must come after including all the getObserverInfo_functions.
#ifdef __cplusplus
static void getObserverDefineMap(ProblemInfo pi, bool clSinglPrecision, std::map<std::string, ObserverInfo> &observerDefineMap, std::vector<std::string> &availableObserverNames) 
{
std::map<std::string, ObserverInfo> newMap;
newMap["basic"]=getObserverInfo_basic(pi, clSinglPrecision);

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
