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

//TODO: expose different observerParams for each observer (provide relevant values only)
//TODO: support using aux vars as event/feature var in observers.
//TODO: concept of "solution buffer" or "solver state" data structure could simplify observer coding. Update it in the ode driver, pass to observer functions

#include "realtype.cl"

#ifdef __cplusplus
template <typename realtype>
#endif
struct ObserverParams
{
    unsigned int eVarIx; //variable for event detection
    unsigned int fVarIx; //variable for features

    unsigned int maxEventCount; //time loop limiter
    unsigned int maxEventTimestamps; //max number of event timestamps to store
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
struct ObserverInfo 
{
	std::string define;
	size_t observerDataSizeFloat;
	size_t observerDataSizeDouble;
	std::vector<std::string> featureNames;
};

#endif

// Design criteria
// - use online algorithms for features such as mean, median, etc.
// - minimize storage of state
// - 

// TODO:
// - return sequence data: list of features per event
// - hashing for distinct value counting [binning?] - not same as above
// - use minimal solution buffer size for task at hand. local extrema: t[3], x[3], dx[2]?
// - separate diagnostic features [toggle on/off independent of observer?]
// - time (e.g., durations) vs state-space features (e.g., amplitudes, means) 
// -- supply list of vars to track for state-space features [similarly for trajectory storage]
// - composable observers? toggle on only what you want 
// - how to handle domain-specific use cases? eg. AHP

// TODO: combine observers with same logic but different event functions into one
// - threshold, nhood (1/2) --> toggle with a define?

////////////////////////////////////////////////
// one-pass detectors
////////////////////////////////////////////////

//basic detectors: no events. Measure extent of state space explored, max/min/mean x and aux, max/min dx
#include "observers/observer_basic.clh" //one variable, specified by fVarIx
#include "observers/observer_basic_allVar.clh"

// convex hull of trajectory in state-space?

// Local maximum detector
// - could do local extremum, toggle whether event is on max vs min?
#include "observers/observer_local_maximum.clh"

// Threshold-based event detection with thresholds defined in state-space coordinates
// #include "observers/observer_threshold_1.clh" //not implemented

// Poincaré section, specified as a normal vector and offset in state-space coordinates
// #include "observers/observer_poincare_1.clh"

// Event trigger is the return of the trajectory to small neighborhood of a point Xstart in state-space coordinates
// - define a sensible Xstart, found in one pass: e.g., local min of a slow variable 
#include "observers/observer_neighborhood_1.clh"

////////////////////////////////////////////////
// two-pass detectors
////////////////////////////////////////////////
// Run a first pass to establish trajectory properties - e.g., extrema for computing normalized state-space coordinates

// Threshold-based event detection with thresholds defined in normalized state-space coordinates
#include "observers/observer_threshold_2.clh"

// Poincaré section, specified as a normal vector and offset in normalized state-space coordinates
// #include "observers/observer_poincare_2.clh"

// Event trigger is the return of the trajectory to small neighborhood of a point Xstart in normalized state-space coordinates
// - Use a first pass to find a good Xstart (e.g. absolute drop below 0.5*range of slowest variable)
#include "observers/observer_neighborhood_2.clh"



// collect available methods into "name"-ObserverInfo map, for C++ side access. Must come after including all the getObserverInfo_functions.
#ifdef __cplusplus
static void getObserverDefineMap(const ProblemInfo pi,
								 const unsigned int fVarIx,
								 const unsigned int eVarIx,
								 const unsigned int nStoredEvents,
								 std::map<std::string, struct ObserverInfo> &observerDefineMap,
								 std::vector<std::string> &availableObserverNames) {
    std::map<std::string, struct ObserverInfo> newMap;
    newMap["basic"]=getObserverInfo_basic(pi, fVarIx, eVarIx, nStoredEvents);
    newMap["basicall"]=getObserverInfo_basicAll(pi, fVarIx, eVarIx, nStoredEvents);
    newMap["localmax"]=getObserverInfo_localmax(pi, fVarIx, eVarIx, nStoredEvents);
    newMap["nhood1"]=getObserverInfo_nhood1(pi, fVarIx, eVarIx, nStoredEvents);
    newMap["nhood2"]=getObserverInfo_nhood2(pi, fVarIx, eVarIx, nStoredEvents);
    newMap["thresh2"]=getObserverInfo_thresh2(pi, fVarIx, eVarIx, nStoredEvents);

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
