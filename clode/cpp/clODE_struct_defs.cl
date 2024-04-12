#ifndef CLODE_STRUCT_DEFS_H_
#define CLODE_STRUCT_DEFS_H_

//TODO: bounds param (like XPP) to catch numerical instability/blow up in finite time?
//TODO: provide different structs for base vs trajectory solvers (expose relevant members)
#include "realtype.cl"

#ifdef __cplusplus
template <typename realtype>
#endif
struct SolverParams
{
	realtype dt;
	realtype dtmax;
	realtype abstol;
	realtype reltol;
	unsigned int max_steps;
	unsigned int max_store;
	unsigned int nout;
};

#endif //CLODE_STRUCT_DEFS_H_
