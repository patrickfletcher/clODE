// header for specifying realtype

//TODO: what if neither is defined?

#ifndef CL_SOLVER_PRECISION_H_
#define CL_SOLVER_PRECISION_H_

//from Sundials package
#if defined(CLODE_SINGLE_PRECISION)

typedef float realtype;
#define RCONST(x) (x##f)
#define BIG_REAL FLT_MAX
#define SMALL_REAL FLT_MIN
#define UNIT_ROUNDOFF FLT_EPSILON

#elif defined(CLODE_DOUBLE_PRECISION)

#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "Double precision floating point not supported by OpenCL implementation."
#endif

typedef double realtype;
#define RCONST(x) (x)
#define BIG_REAL DBL_MAX
#define SMALL_REAL DBL_MIN
#define UNIT_ROUNDOFF DBL_EPSILON

#endif

#endif //CL_SOLVER_PRECISION_H_
