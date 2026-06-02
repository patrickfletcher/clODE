#ifndef OBSERVERS_H_
#define OBSERVERS_H_

/* "Observer" measures features of the ODE solution as it is being integrated
 * The observer consists of a persistent state structure and several functions:
 * 
 * - initializeObserverState: seed persistent state and retained event/output buffers
 * - warmupObserverState: for two-pass detectors, collect only warmup-pass geometry
 * - updateObserverState: advance accepted-step history buffers and continuous per-step reducers
 * - initializeEventDetector: finalize live-pass detector geometry using warmup data and the rewound initial sample
 * - eventFunction: decide whether the observer's public event fired on the updated history, optionally refining event geometry within the step
 * - computeEventFeatures: when an event is detected, update retained outputs and eventwise aggregates such as counts or periods
 * - finalizeFeatures: post-integration cleanup and write to global feature array
 */

//TODO: expose different observer runtime settings for each observer (provide relevant values only)
//TODO: support using aux vars as event/feature var in observers.
//TODO: concept of "solution buffer" or "solver state" structure could simplify observer coding. Update it in the ode driver, pass to observer functions

#include "clODE_struct_defs.cl"

// Shared accepted-step history helpers for observer-local solution buffers.
// Keep scalar K=2/K=3 helpers unrolled for the per-step hot path.
static inline void advanceAcceptedStepHistory2(realtype *history, const realtype newestSample) {
	history[0] = history[1];
	history[1] = newestSample;
}

// Update a packed per-variable K=2 history laid out as [var0_k0, var0_k1, var1_k0, ...].
static inline void advanceAcceptedStepHistory2ByVariable(realtype *packedHistory, const realtype latestValues[], const int valueCount) {
	for (int j = 0; j < valueCount; ++j) {
		realtype *history = &packedHistory[j * 2];
		advanceAcceptedStepHistory2(history, latestValues[j]);
	}
}

static inline void advanceAcceptedStepHistory3(realtype *history, const realtype newestSample) {
	history[0] = history[1];
	history[1] = history[2];
	history[2] = newestSample;
}

// Update a packed per-variable K=3 history laid out as [var0_k0, var0_k1, var0_k2, var1_k0, ...].
// The loop-based form keeps call sites compact when advancing N_VAR-sized buffers.
static inline void advanceAcceptedStepHistory3ByVariable(realtype *packedHistory, const realtype latestValues[], const int valueCount) {
	for (int j = 0; j < valueCount; ++j) {
		realtype *history = &packedHistory[j * 3];
		advanceAcceptedStepHistory3(history, latestValues[j]);
	}
}

static inline realtype observerElapsedValue(
	const realtype elapsedAtChunkStart,
	const realtype elapsedAtChunkStartCorrection,
	const realtype solveElapsedTotal
) {
	return compensatedTimeValueAfterStep(
		elapsedAtChunkStart,
		elapsedAtChunkStartCorrection,
		solveElapsedTotal
	);
}

static inline void advanceObserverElapsedChunk(
	realtype *elapsedAtChunkStart,
	realtype *elapsedAtChunkStartCorrection,
	const realtype solveElapsedTotal
) {
	compensatedTimeAdd(
		elapsedAtChunkStart,
		elapsedAtChunkStartCorrection,
		solveElapsedTotal
	);
}

// should replicate:
// for (int j = 0; j < N_VAR; ++j) {
// 	od->xbuffer[j * 3 + 0] = od->xbuffer[j * 3 + 1];
// 	od->xbuffer[j * 3 + 1] = od->xbuffer[j * 3 + 2];
// 	od->xbuffer[j * 3 + 2] = xi[j];}


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
#include "observers/observer_threshold_crossing.clh"

// Schmitt-trigger event detection with absolute thresholds defined in state-space coordinates
#include "observers/observer_schmitt_trigger.clh"

// Poincaré section, specified as a normal vector and offset in state-space coordinates
// #include "observers/observer_poincare_1.clh"

// Event trigger is the return of the trajectory to small neighborhood of a point Xstart in state-space coordinates
// - define a sensible Xstart, found in one pass: e.g., local min of a slow variable 
// #include "observers/observer_neighborhood_1.clh"

////////////////////////////////////////////////
// two-pass detectors
////////////////////////////////////////////////
// Run a first pass to establish trajectory properties - e.g., extrema for computing normalized state-space coordinates

// Threshold-based event detection with one warmup-derived normalized boundary
#include "observers/observer_normalized_threshold_crossing.clh"

// Schmitt-trigger event detection with warmup-derived normalized thresholds
#include "observers/observer_normalized_schmitt_trigger.clh"

// Poincaré section, specified as a normal vector and offset in normalized state-space coordinates
// #include "observers/observer_poincare_2.clh"

// Event trigger is the return of the trajectory to a small normalized neighborhood of a warmup-derived anchor point
#include "observers/observer_normalized_neighborhood_return.clh"


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
