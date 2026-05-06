//collection of helper functions useful in feature detectors/RHS function computations.

// TODO: extended precision helpers (TwoSum etc) - in realtype.cl?
// TODO: attributes(e.g., align, static inline, ...), #pragma unroll, etc?
// TOOD: compiler flags (fma, mad, ...?)
// - e.g. DEFINE to swap in alternatives? 

// TODO: expand interpolation routines [use slope info to provide better accuracy]
// - DEFINE to swap in method alternatives [none, linear, quad, etc]


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

//TODO 2sum for accurate summation (e.g. t+=dt...)?
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm

// Neumaier’s algorithm for summation - avoid loss of precision in accumulators
// - requires adding the correction at the end, once
// static inline void twoSum(realtype *sum, realtype *correction, realtype newValue)
// {
// 	realtype total;
// 	volatile realtype diff; //is needed?
// 	t = *sum + newValue;
// 	if (fabs(*sum) >= fabs(newValue))
// 		diff = (*sum - t);
// 		*correction += diff + newValue;
// 	else 
// 		diff = (newValue - t);
// 		*correction += diff + *sum;
// 	*sum = t;
// }

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

//estimate yi at specified ti, using linear interpolation of two values
static inline realtype linearInterp(realtype t0, realtype t1, realtype y0, realtype y1, realtype ti) {
	realtype yi = y0 + (ti - t0) * (y1 - y0) / (t1 - t0);
	return yi;
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

//compute vertex of a quadratic interpolant of three values
// - store result in tv, yv
static inline void quadraticInterpVertex(realtype t[], realtype y[], realtype *tv, realtype *yv) {
	realtype b0, b1, b2;

	b0 = y[0];
	b1 = (y[1] - b0) / (t[1] - t[0]);
	b2 = (y[2] - b0 - b1 * (t[2] - t[0])) / ((t[2] - t[0]) * (t[2] - t[1]));

	*tv = -(b1 - b2 * (t[0] + t[1])) / (RCONST(2.0) * b2);
	*yv = b0 + b1 * (*tv - t[0]) + b2 * (*tv - t[0]) * (*tv - t[1]);
}

//~ static inline realtype cubicInterp(realtype t[], realtype y[], realtype dy[], realtype ti) {
//~ return yi;
//~ }

#endif //CL_UTILITIES_H_
