#ifndef CLODE_STRUCT_DEFS_H_
#define CLODE_STRUCT_DEFS_H_

//TODO: bounds param (like XPP) to catch numerical instability/blow up in finite time?
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
	int max_steps;
	int max_store;
	int nout;
};

#endif //CLODE_STRUCT_DEFS_H_
