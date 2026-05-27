#ifndef OBSERVERS_H_
#define OBSERVERS_H_

/* "Observer" measures features of the ODE solution as it is being integrated
 * The observer consists of a persistent state structure and several functions:
 * 
 * - initializeObserverState: set up the observer state to sensible values
 * - warmupObserverState: for two-pass event detectors - restricted data collection about trajectory during a first pass ODE solve
 * - updateObserverState: per-timestep update of observer state
 * - initializeEventDetector: set any values needed to do selected type of event detection (possibly using warmup data)
 * - eventFunction: check for an event. Optionally refine location of event within timestep. Compute event-based quantities
 * - computeEventFeatures: when event is detected, compute desired per-event features
 * - finalizeFeatures: post-integration cleanup and write to global feature array
 */

//TODO: expose different observer runtime settings for each observer (provide relevant values only)
//TODO: support using aux vars as event/feature var in observers.
//TODO: concept of "solution buffer" or "solver state" structure could simplify observer coding. Update it in the ode driver, pass to observer functions

#include "clODE_struct_defs.cl"


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

//basic detectors: no events. Measure selected max/min/mean x and aux, max/min/mean dx
#include "observers/observer_summary.clh"

// convex hull of trajectory in state-space?

// Local maximum detector
// - could do local extremum, toggle whether event is on max vs min?
#include "observers/observer_local_maximum.clh"

// Threshold-based event detection with thresholds defined in state-space coordinates
#include "observers/observer_threshold_1.clh"

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
