//TODO: Force all steppers to follow FSAL? allows "free" interpolation in any timestep (important for event detection, trajectory output at specified times).
//      if purifying ti (fixed steppers) must get new dxi in main loop.
//TODO: Adaptive stepper seems to fail to keep stepsize small enough to prevent blowups. Something to do with abstol??
//TODO: stepper wishlist: Multi-step method support,  implicit methods
//TODO: does fma(a,b,c)=a*b+c for improved accuracy? vs mad? Does compiler do it anyway (if so write normal version for readability)?
//TODO: fixed stepsize: if always purifying ti in main loop (t0+step*dt), don't update ti here?
//TODO: test speed increases for "minimizing registers" vs reuse of computations (e.g. h2=*dt*0.5)

#ifndef STEPPERS_H_
#define STEPPERS_H_

//available methods "manifest" for C++ use
#ifdef __cplusplus 

#include <string>
#include <map>
#include <vector>

//store name/define pairs of available steppers as a map for access in C++.  Static when defining function in header
static void getStepperDefineMap(std::map<std::string, std::string> &stepperDefineMap, std::vector<std::string> &availableStepperNames) 
{
std::map<std::string, std::string> newMap;
newMap["euler"]="EXPLICIT_EULER";
newMap["heun"]="EXPLICIT_HEUN";
newMap["rk4"]="EXPLICIT_RK4";
newMap["bs23"]="EXPLICIT_BS23";
newMap["dopri5"]="EXPLICIT_DOPRI5";
newMap["seuler"]="STOCHASTIC_EULER";

//export vector of names for access in C++
std::vector<std::string> newNames;
for (auto const& element : newMap)
    newNames.push_back(element.first);

stepperDefineMap=newMap;
availableStepperNames=newNames;
}

#else

#include "realtype.cl"

//forward declaration of the RHS function
void getRHS(const realtype t, const realtype x_[], const realtype p_[], realtype dx_[], realtype aux_[], const realtype w_[]);

// FIXED STEPSIZE EXPLICIT METHODS
// The fixed steppers use stepcount to purify the T values (eliminates roundoff)

#ifdef STOCHASTIC_EULER
#define STOCHASTIC_STEPPER
#include "steppers/fixed_explicit_Euler.clh"
#endif

#ifdef EXPLICIT_EULER
#include "steppers/fixed_explicit_Euler.clh"
#endif

#ifdef EXPLICIT_HEUN
#include "steppers/fixed_explicit_Trapezoidal.clh"
#endif

//~ #ifdef MIDPOINT
//~ #endif

#ifdef EXPLICIT_RK4
#include "steppers/fixed_explicit_RK4.clh"
#endif

//~ #ifdef RK higher order?
//~ #endif

//~ #ifdef SSPRK3
//~ #endif

#ifdef FIXED_STEPSIZE_EXPLICIT
#include "steppers/fixed_explicit_step.clh"
#endif

// TODO: FIXED STEP MULTI STEP METHODS
// need to fill the first few steps with a single stepper, so if MULTI_STEP is defined also include Euler or RK4 or something
// #define MULTI_STEP


// ADAPTIVE STEPSIZE EXPLICIT METHODS
#ifdef HEUN_EULER
#include "steppers/adaptive_he12.clh"
#endif

#ifdef EXPLICIT_BS23
#include "steppers/adaptive_bs23.clh"
#endif

//~ #ifdef CASH_KARP
//~ #endif

#ifdef EXPLICIT_DOPRI5
#include "steppers/adaptive_dp45.clh"
#endif

//~ #ifdef RKF87
//~ #endif

// low-register use adaptive stepsize solvers, from rktides
//~ #ifdef RK435_KENNEDY
//~ #endif

//~ #ifdef RK549_KENNEDY
//~ #endif

#ifdef ADAPTIVE_STEPSIZE_EXPLICIT
#include "steppers/adaptive_explicit_step.clh"
#endif



//TODO: implicit fixed/adaptive steppers




#endif //__cplusplus
#endif //STEPPERS_H_