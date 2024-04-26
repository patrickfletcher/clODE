// header for specifying realtype: precision and related items

// TODO: extended precision helpers (TwoSum)

#ifndef CL_SOLVER_PRECISION_H_
#define CL_SOLVER_PRECISION_H_

//similar to Sundials package   -------------
#if defined(CLODE_SINGLE_PRECISION)

typedef float realtype;
typedef float2 realtype2;
typedef float4 realtype4;
typedef float8 realtype8;
typedef float16 realtype16;
#define RCONST(x) (x##f)
#define BIG_REAL FLT_MAX
#define SMALL_REAL FLT_MIN
#define UNIT_ROUNDOFF FLT_EPSILON
#define ZERO 0.0f
#define ONE 1.0f

#elif defined(CLODE_DOUBLE_PRECISION) //CLODE_SINGLE_PRECISION  -------------

#ifdef cl_khr_fp64
#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#elif defined(cl_amd_fp64)
#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#else
#error "Double precision floating point not supported by OpenCL implementation."
#endif //cl_khr_fp64 -------------

typedef double realtype;
typedef double2 realtype2;
typedef double4 realtype4;
typedef double8 realtype8;
typedef double16 realtype16;
#define RCONST(x) (x)
#define BIG_REAL DBL_MAX
#define SMALL_REAL DBL_MIN
#define UNIT_ROUNDOFF DBL_EPSILON
#define ZERO 0.0
#define ONE 1.0

#endif //CLODE_DOUBLE_PRECISION  -------------

#endif //CL_SOLVER_PRECISION_H_
