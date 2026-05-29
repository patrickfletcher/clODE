//collection of helper functions useful in feature detectors/RHS function computations.

// TODO: extended precision helpers (TwoSum etc) - in realtype.cl?
// TODO: attributes(e.g., align, static inline, ...), #pragma unroll, etc?
// TOOD: compiler flags (fma, mad, ...?)
// - e.g. DEFINE to swap in alternatives? 

// TODO: expand interpolation routines [use slope info to provide better accuracy]
// - DEFINE to swap in method alternatives [none, linear, quad, etc]

// NOTES:
// - work-item-local arrays and struct members are __private in OpenCL C.
// - array helpers here are currently intended for caller-owned private state or
//   scratch storage unless documented otherwise.
// - if a shared helper needs buffer-backed data, prefer an explicit named
//   address-space signature (__global/__local/__constant) instead of relying on
//   __generic or implementation-specific defaults.
// - can pass in N_VAR as a parameter if needed, see for example norm_inf.

#ifndef CL_UTILITIES_H_
#define CL_UTILITIES_H_

#include "realtype.cl"

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define heaviside(x) ((x) >= ZERO ? ONE : ZERO) 

//1-norm
static inline realtype norm_1(realtype x[], int N) {
	realtype result = ZERO;
	for (int k = 0; k < N; k++)
		result += fabs(x[k]);

	return result;
}

//2-norm
static inline realtype norm_2(realtype x[], int N) {
	realtype result = ZERO;
	for (int k = 0; k < N; k++)
		result += x[k] * x[k];

	return sqrt(result);
}

//Inf-norm
static inline realtype norm_inf(realtype x[], int N) {
	realtype result = ZERO;
	for (int k = 0; k < N; k++)
		result = fmax(fabs(x[k]), result);

	return result;
}

//Maximum of vector, returns both max and index of max
static inline void maxOfArray(realtype inArray[], int N, realtype *maxVal, int *index) {
	*maxVal = -BIG_REAL;
	*index = 0;
	for (int k = 0; k < N; k++)
	{
		if (inArray[k] > *maxVal) //return first occurrence (>= returns last)
		{
			*maxVal = inArray[k];
			*index = k;
		}
	}
}

// returns maximum value in a 1D array
static inline realtype array_max(realtype inArray[], int N) {
	realtype maxVal = -BIG_REAL;
	for (int k = 0; k < N; k++)
	{
		if (inArray[k] > maxVal) //return first occurrence (>= returns last)
			maxVal = inArray[k];
	}
	return maxVal;
}

// returns index of maximum value in a 1D array
static inline int array_argmax(realtype inArray[], int N) {
	realtype maxVal = -BIG_REAL;
	int index = 0;
	for (int k = 0; k < N; k++)
	{
		if (inArray[k] > maxVal) //return first occurrence (>= returns last)
		{
			maxVal = inArray[k];
			index = k;
		}
	}
	return index;
}


//Minimum of vector, returns both min and index of min
static inline void minOfArray(realtype inArray[], int N, realtype *minVal, int *index) {
	*minVal = BIG_REAL;
	*index = 0;
	for (int k = 0; k < N; k++)
	{
		if (inArray[k] < *minVal) //return first occurrence (<= returns last)
		{
			*minVal = inArray[k];
			*index = k;
		}
	}
}

// returns minimum value in a 1D array
static inline realtype array_min(realtype inArray[], int N) {
	realtype minVal = BIG_REAL;
	for (int k = 0; k < N; k++)
	{
		if (inArray[k] < minVal) //return first occurrence (<= returns last)
			minVal = inArray[k];
	}
	return minVal;
}

// returns index of minimum value in a 1D array
static inline int array_argmin(realtype inArray[], int N) {
	realtype minVal = BIG_REAL;
	int index = 0;
	for (int k = 0; k < N; k++)
	{
		if (inArray[k] < minVal) //return first occurrence (<= returns last)
		{
			minVal = inArray[k];
			index = k;
		}
	}
	return index;
}


/* 
Online algorithms for feature detection
Goal: compute features using minimal storage as the solution is numerically approximated
- running means (should handle non-uniform sampling for adaptive time-stepping)
-- could implement an online trapezoidal method for running mean?
- running percentiles - P2 algorithm

other notes
- may be of interest to use larger solution buffer (>3) along with some notion of "forgetting" old states?
- check, e.g., online/real-time speech/signal processing literature
*/

// compensated summation

// Neumaier’s algorithm for summation - avoid loss of precision in accumulators.
// The caller owns both the running sum and the correction term so kernels can
// adopt compensation selectively without introducing a larger state object.
static inline void compensatedSumAdd(realtype *sum, realtype *correction, realtype newValue) {
	realtype total = *sum + newValue;
	if (fabs(*sum) >= fabs(newValue))
		*correction += (*sum - total) + newValue;
	else
		*correction += (newValue - total) + *sum;
	*sum = total;
}

static inline realtype compensatedSumValue(realtype sum, realtype correction) {
	return sum + correction;
}

static inline void compensatedIntegrateConstant(
	realtype *integral,
	realtype *correction,
	realtype dt,
	realtype value
) {
	if (dt <= ZERO)
		return;
	compensatedSumAdd(integral, correction, dt * value);
}

static inline realtype meanFromCompensatedIntegral(
	realtype integral,
	realtype correction,
	realtype total_delta
) {
	if (total_delta <= ZERO)
		return ZERO;
	return compensatedSumValue(integral, correction) / total_delta;
}

// Kahan-style time accumulation for solver-relative elapsed-time bookkeeping.
// The correction term stores the low-order bits that were lost by the running
// sum, so reconstructing the compensated value subtracts that term back out.
// This mutates the elapsed time and correction in place - performs the true updates to timeValue and timeCorrection
// Reference: https://en.wikipedia.org/wiki/Kahan_summation_algorithm
static inline void compensatedTimeAdd(realtype *timeValue, realtype *timeCorrection, realtype dt) {
	if (dt == ZERO)
		return;
	realtype y = dt - *timeCorrection;
	realtype total = *timeValue + y;
	*timeCorrection = (total - *timeValue) - y;
	*timeValue = total;
}

// reconstruct the compensated time value from the running sum and correction term - does not mutate either term
static inline realtype compensatedTimeValue(realtype timeValue, realtype timeCorrection) {
	return timeValue - timeCorrection;
}

// compute the compensated time value after adding a step, without mutating the running sum or correction term
static inline realtype compensatedTimeValueAfterStep(
	realtype timeValue,
	realtype timeCorrection,
	realtype dt
) {
	compensatedTimeAdd(&timeValue, &timeCorrection, dt);
	return compensatedTimeValue(timeValue, timeCorrection);
}

// compute the compensated time from the origin, given the current elapsed time and correction - does not mutate either term
static inline realtype compensatedTimeFromOrigin(
	realtype origin,
	realtype elapsed,
	realtype elapsedCorrection
) {
	return origin + compensatedTimeValue(elapsed, elapsedCorrection);
}

// this is read-only; it doesn't update the elapsed time or correction, so it can be used to compute stage times without mutating the state of the time accumulator
static inline realtype compensatedTimeFromOriginAfterStep(
	realtype origin,
	realtype elapsed,
	realtype elapsedCorrection,
	realtype dt
) {
	return origin + compensatedTimeValueAfterStep(elapsed, elapsedCorrection, dt);
}

// TODO: evaluate incremental versions (below) vs running sum (two-sum) then a single division at the end. Need to do so for variance already anyway

//Compute a running mean of a function at possibly non-uniform sample points
// - mean should be initialized to zero externally (first step: dt=total_delta --> mean=newValue)
static inline realtype runningMeanTime(realtype mean, realtype newValue, realtype dt, realtype total_delta) {
	return mean + (newValue - mean) * dt/total_delta;
}

//Compute a running mean for a set of numbers
static inline void runningMean(realtype *mean, realtype newValue, unsigned int eventCount) {
	if (eventCount == 1) //initialize the mean to the first value
		*mean = newValue;
	else if (eventCount > 1) //compute the current value of the running mean
		*mean += (newValue - *mean) / (realtype)eventCount;
}

//Compute a running mean and variance for a set of numbers
// https://www.johndcook.com/blog/standard_deviation/
// NOTE: once the variance value is desired, it must be divided by the final event count!
static inline void runningMeanVar(realtype *mean, realtype *variance, realtype newValue, unsigned int eventCount) {
	if (eventCount == 1)
	{ //initialize the mean to the first value, variance to zero
		*mean = newValue;
		*variance = ZERO;
	}
	else if (eventCount > 1)
	{ //compute the current value of the running mean and variance
		realtype tmp = *mean;
		*mean = tmp + (newValue - tmp) / (realtype)eventCount;
		*variance = *variance + (newValue - tmp) * (newValue - *mean);
	}
}


// Interpolation routines
// TODO: try quadratic interp using 3 points for threshold crossings too.

//estimate yi at specified ti, using linear interpolation of two values
static inline realtype linearInterp(realtype t0, realtype t1, realtype y0, realtype y1, realtype ti) {
	realtype yi = y0 + (ti - t0) * (y1 - y0) / (t1 - t0);
	return yi;
}

// Estimate the time at which y reaches yi between two samples.
static inline realtype linearInterpTimeOfValue(realtype t0, realtype t1, realtype y0, realtype y1, realtype yi) {
	if (t1 == t0 || y1 == y0)
		return t1;
	return t0 + (yi - y0) * (t1 - t0) / (y1 - y0);
}

//estimate yi at specified ti, using linear interpolation between the first or second pair of values, given three values 
// - the solution buffer in clode keeps t/y values of the most recent 3 time steps
static inline realtype linearInterpArray(realtype t[], realtype y[], realtype ti) {
	realtype yi;
	if (ti < t[1])
		yi = y[0] + (ti - t[0]) * (y[1] - y[0]) / (t[1] - t[0]);
	else
		yi = y[1] + (ti - t[1]) * (y[2] - y[1]) / (t[2] - t[1]);

	return yi;
}

//estimate yi at specified ti, using quadratic interpolant of three values
static inline realtype quadraticInterp(realtype t[], realtype y[], realtype ti) {
	realtype b0, b1, b2, yi;

	b0 = y[0];
	b1 = (y[1] - b0) / (t[1] - t[0]);
	b2 = (y[2] - b0 - b1 * (t[2] - t[0])) / ((t[2] - t[0]) * (t[2] - t[1]));

	yi = b0 + b1 * (ti - t[0]) + b2 * (ti - t[0]) * (ti - t[1]);

	return yi;
}

static inline realtype cubicHermiteValueUnitInterval(
	realtype u,
	realtype y0,
	realtype y1,
	realtype dy0,
	realtype dy1,
	realtype dt
) {
	realtype h00 = RCONST(2.0) * u * u * u - RCONST(3.0) * u * u + ONE;
	realtype h10 = u * u * u - RCONST(2.0) * u * u + u;
	realtype h01 = -RCONST(2.0) * u * u * u + RCONST(3.0) * u * u;
	realtype h11 = u * u * u - u * u;
	return h00 * y0 + h10 * dt * dy0 + h01 * y1 + h11 * dt * dy1;
}

static inline realtype cubicHermiteDerivativeUnitInterval(
	realtype u,
	realtype y0,
	realtype y1,
	realtype dy0,
	realtype dy1,
	realtype dt
) {
	realtype dh00 = RCONST(6.0) * u * u - RCONST(6.0) * u;
	realtype dh10 = RCONST(3.0) * u * u - RCONST(4.0) * u + ONE;
	realtype dh01 = -RCONST(6.0) * u * u + RCONST(6.0) * u;
	realtype dh11 = RCONST(3.0) * u * u - RCONST(2.0) * u;
	return dh00 * y0 + dh10 * dt * dy0 + dh01 * y1 + dh11 * dt * dy1;
}

// Prototype only: use endpoint values and slopes to refine a threshold-crossing
// time inside one timestep. Falls back to linear inversion if Newton leaves the
// unit interval or the cubic derivative becomes too small.
// Reference: https://en.wikipedia.org/wiki/Cubic_Hermite_spline
static inline realtype cubicHermiteInterpTimeOfValue(
	realtype t0,
	realtype t1,
	realtype y0,
	realtype y1,
	realtype dy0,
	realtype dy1,
	realtype yi
) {
	if (t1 == t0 || y1 == y0)
		return t1;

	realtype dt = t1 - t0;
	realtype linearGuess = linearInterpTimeOfValue(t0, t1, y0, y1, yi);
	realtype u = clamp((linearGuess - t0) / dt, ZERO, ONE);

	for (int iter = 0; iter < 4; ++iter) {
		realtype value = cubicHermiteValueUnitInterval(u, y0, y1, dy0, dy1, dt);
		realtype deriv = cubicHermiteDerivativeUnitInterval(u, y0, y1, dy0, dy1, dt);
		if (fabs(deriv) <= UNIT_ROUNDOFF)
			return linearGuess;
		realtype candidate = u - (value - yi) / deriv;
		if (candidate <= ZERO || candidate >= ONE)
			return linearGuess;
		u = candidate;
	}

	return t0 + u * dt;
}

//compute vertex of a quadratic interpolant of three values
// TODO: consider detection of extrema using 3 x values is better than dx sign change?
// - store result in tv, yv
static inline void quadraticInterpVertex(realtype t[], realtype y[], realtype *tv, realtype *yv) {
	realtype b0, b1, b2;

	b0 = y[0];
	b1 = (y[1] - b0) / (t[1] - t[0]);
	b2 = (y[2] - b0 - b1 * (t[2] - t[0])) / ((t[2] - t[0]) * (t[2] - t[1]));

	*tv = -(b1 - b2 * (t[0] + t[1])) / (RCONST(2.0) * b2);
	*yv = b0 + b1 * (*tv - t[0]) + b2 * (*tv - t[0]) * (*tv - t[1]);
}

static inline bool quadraticInterpVertexBounded(realtype t[], realtype y[], realtype *tv, realtype *yv) {
	realtype dt10 = t[1] - t[0];
	realtype dt20 = t[2] - t[0];
	realtype dt21 = t[2] - t[1];
	if (dt10 == ZERO || dt20 == ZERO || dt21 == ZERO)
		return false;

	realtype b0 = y[0];
	realtype b1 = (y[1] - b0) / dt10;
	realtype b2 = (y[2] - b0 - b1 * dt20) / (dt20 * dt21);
	if (fabs(b2) <= UNIT_ROUNDOFF)
		return false;

	realtype candidateT = -(b1 - b2 * (t[0] + t[1])) / (RCONST(2.0) * b2);
	if (candidateT < t[0] || candidateT > t[2])
		return false;

	realtype candidateY = b0 + b1 * (candidateT - t[0]) + b2 * (candidateT - t[0]) * (candidateT - t[1]);
	*tv = candidateT;
	*yv = candidateY;
	return true;
}

static inline void localMaximumFromThreeSamples(realtype t[], realtype y[], realtype *tMax, realtype *yMax) {
	int index = array_argmax(y, 3);
	*tMax = t[index];
	*yMax = y[index];

	realtype candidateT;
	realtype candidateY;
	if (quadraticInterpVertexBounded(t, y, &candidateT, &candidateY) && candidateY >= *yMax) {
		*tMax = candidateT;
		*yMax = candidateY;
	}
}

static inline void localMinimumFromThreeSamples(realtype t[], realtype y[], realtype *tMin, realtype *yMin) {
	int index = array_argmin(y, 3);
	*tMin = t[index];
	*yMin = y[index];

	realtype candidateT;
	realtype candidateY;
	if (quadraticInterpVertexBounded(t, y, &candidateT, &candidateY) && candidateY <= *yMin) {
		*tMin = candidateT;
		*yMin = candidateY;
	}
}

//~ static inline realtype cubicInterp(realtype t[], realtype y[], realtype dy[], realtype ti) {
//~ return yi;
//~ }

#endif //CL_UTILITIES_H_
