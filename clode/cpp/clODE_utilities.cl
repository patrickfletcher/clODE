//collection of helper functions useful in feature detectors/RHS function computations.

#ifndef CL_UTILITIES_H_
#define CL_UTILITIES_H_

#include "realtype.cl"

//computation done by preprocessor:
#define MIN(a, b) ((a) < (b) ? (a) : (b)) //TODO: use builtin fmin/fmax (just make sure to not mix arg type)
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define heav(x) ((x) > RCONST(0.0) ? RCONST(1.0) : RCONST(0.0)) //TODO: capitalize? Or use OpenCL step (just make sure to not mix arg type)

//TODO 2sum for accurate summation (e.g. t+=dt...)?
// https://en.wikipedia.org/wiki/Kahan_summation_algorithm

//1-norm
inline realtype norm_1(realtype x[], int N)
{
	realtype result = RCONST(0.0);
	for (int k = 0; k < N; k++)
		result += fabs(x[k]);

	return result;
}

//2-norm
inline realtype norm_2(realtype x[], int N)
{
	realtype result = RCONST(0.0);
	for (int k = 0; k < N; k++)
		result += x[k] * x[k];

	return sqrt(result);
}

//Inf-norm
inline realtype norm_inf(realtype x[], int N)
{
	realtype result = RCONST(0.0);
	for (int k = 0; k < N; k++)
		result = fmax(fabs(x[k]), result);

	return result;
}

//Maximum of vector, returns both max and index of max
inline void maxOfArray(realtype inArray[], int N, realtype *maxVal, int *index)
{
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

//Minimum of vector, returns both min and index of min
inline void minOfArray(realtype inArray[], int N, realtype *minVal, int *index)
{
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

//Compute a running mean
inline void runningMean(realtype *mean, realtype thisValue, unsigned int eventCount)
{
	if (eventCount == 0)
	{ //do nothing
	}
	else if (eventCount == 1)
	{ //initialize the mean to the first value
		*mean = thisValue;
	}
	else
	{ //compute the current value of the running mean
		*mean = *mean + (thisValue - *mean) / (realtype)eventCount;
	}
}

//Compute a running mean and variance. 
// https://www.johndcook.com/blog/standard_deviation/
// NOTE: once the variance value is desired, it must be divided by the final event count!
inline void runningMeanVar(realtype *mean, realtype *variance, realtype thisValue, unsigned int eventCount)
{
	if (eventCount == 0)
	{ //do nothing
	}
	else if (eventCount == 1)
	{ //initialize the mean to the first value, variance to zero
		*mean = thisValue;
		*variance = RCONST(0.0);
	}
	else
	{ //compute the current value of the running mean and variance
		realtype tmp = *mean;
		*mean = tmp + (thisValue - tmp) / (realtype)eventCount;
		*variance = *variance + (thisValue - tmp) * (thisValue - *mean);
	}
}

//estimate yi at specified ti, using linear interpolation of two values
inline realtype linearInterp(realtype t0, realtype t1, realtype y0, realtype y1, realtype ti)
{
	realtype yi = y0 + (ti - t0) * (y1 - y0) / (t1 - t0);
	return yi;
}

//estimate yi at specified ti, using linear interpolation between the first or second pair of values, given three values 
// - the solution buffer in clode keeps t/y values of the most recent 3 time steps
inline realtype linearInterpArray(realtype t[], realtype y[], realtype ti)
{
	realtype yi;
	if (ti < t[1])
		yi = y[0] + (ti - t[0]) * (y[1] - y[0]) / (t[1] - t[0]);
	else
		yi = y[1] + (ti - t[1]) * (y[2] - y[1]) / (t[2] - t[1]);

	return yi;
}

//estimate yi at specified ti, using quadratic interpolant of three values
inline realtype quadraticInterp(realtype t[], realtype y[], realtype ti)
{
	realtype b0, b1, b2, yi;

	b0 = y[0];
	b1 = (y[1] - b0) / (t[1] - t[0]);
	b2 = (y[2] - b0 - b1 * (t[2] - t[0])) / ((t[2] - t[0]) * (t[2] - t[1]));

	yi = b0 + b1 * (ti - t[0]) + b2 * (ti - t[0]) * (ti - t[1]);

	return yi;
}

//compute vertex of a quadratic interpolant of three values
// - store result in tv, yv
inline void quadraticInterpVertex(realtype t[], realtype y[], realtype *tv, realtype *yv)
{
	realtype b0, b1, b2;

	b0 = y[0];
	b1 = (y[1] - b0) / (t[1] - t[0]);
	b2 = (y[2] - b0 - b1 * (t[2] - t[0])) / ((t[2] - t[0]) * (t[2] - t[1]));

	*tv = -(b1 - b2 * (t[0] + t[1])) / (RCONST(2.0) * b2);
	*yv = b0 + b1 * (*tv - t[0]) + b2 * (*tv - t[0]) * (*tv - t[1]);
}

//TODO: use slopes!
// quadratic given slopes, cubic given slopes, ...
//~ inline realtype cubicInterp(realtype t[], realtype y[], realtype dy[], realtype ti) {
//~ return yi;
//~ }

#endif //CL_UTILITIES_H_
