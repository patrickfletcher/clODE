//TODO: Force all steppers to follow FSAL? allows "free" interpolation in any timestep (important for event detection, trajectory output at specified times).
//      if purifying ti (fixed steppers) must get new dxi in main loop.
//TODO: Adaptive stepper seems to fail to keep stepsize small enough to prevent blowups. Something to do with abstol??
//TODO: stepper wishlist: Multi-step method support,  implicit methods
//TODO: does fma(a,b,c)=a*b+c for improved accuracy? vs mad? Does compiler do it anyway (if so write normal version for readability)?
//TODO: fixed stepsize: if always purifying ti in main loop (t0+step*dt), don't update ti here?
//TODO: test speed increases for "minimizing registers" vs reuse of computations (e.g. h2=*dt*0.5)

#ifndef STEPPERS_H_
#define STEPPERS_H_

#include "clODE_struct_defs.cl" //for SolverParams struct definition
#include "realtype.cl"

//forward declaration of the RHS function
void getRHS(realtype t, realtype x_[], realtype p_[], realtype dx_[], realtype aux_[], realtype w_[]);

//store name/define pairs of available steppers as a map for access in C++
#ifdef __cplusplus
#include <unordered_map>
std::unordered_map<std::string,std::string> stepperDefineMap;
stepperDefineMap["explicitEuler"]="EXPLICIT_EULER";
stepperDefineMap["explicitTrapezoidal"]="EXPLICIT_TRAPEZOIDAL";
#endif

// FIXED STEPSIZE EXPLICIT METHODS
/* The fixed steppers use stepcount to purify the T values (eliminates roundoff)
 * This means the returned "k1" is dx at the previous step; need one extra DX computation for the final step
 */
#ifdef FORWARD_EULER
#include "steppers/fixed_explicitEuler.clh"
#endif

#ifdef HEUN
#include "steppers/fixed_explicitTrapezoidal.clh"
#endif

//~ #ifdef MIDPOINT
//~ #endif

#ifdef RUNGE_KUTTA4
#include "steppers/fixed_explicitRK4.clh"
#endif

//~ #ifdef RK higher order?
//~ #endif

//~ #ifdef SSPRK3
//~ #endif

// TODO: FIXED STEP MULTI STEP METHODS
// need to fill the first few steps with a single stepper, so if MULTI_STEP is defined also include Euler or RK4 or something
// #define MULTI_STEP


// ADAPTIVE STEPSIZE EXPLICIT METHODS
#ifdef HEUN_EULER
#include "steppers/adaptive_he12.clh"
#endif

#ifdef BOGACKI_SHAMPINE23
#include "steppers/adaptive_bs23.clh"
#endif

//~ #ifdef CASH_KARP
//~ #endif

#ifdef DORPRI5
#include "steppers/adaptive_dp45.clh"
#endif

//~ #ifdef RKF87
//~ #endif

// low-register use adaptive stepsize solvers, from rktides
//~ #ifdef RK435_KENNEDY
//~ #endif

//~ #ifdef RK549_KENNEDY
//~ #endif


//TODO: implicit fixed/adaptive steppers


#ifdef ADAPTIVE_STEPSIZE

//Wrapper to handle step-size adaptation.  note: wi should be zeros
int stepper(realtype *ti, realtype xi[], realtype k1[], realtype pars[], __constant struct SolverParams *sp, realtype *dt, __constant realtype *tspan, realtype aux[], realtype wi[])
{
    realtype err[N_VAR];
    realtype normErr, newDt = *dt;
    realtype tNew, newxi[N_VAR], newk1[N_VAR];

    realtype threshold = sp->abstol / sp->reltol;
    realtype expon = RCONST(1.0) / ORDERPLUSONE;
    realtype hmin = RCONST(16.0) * fabs(fabs(nextafter(*ti, RCONST(1.1)*tspan[1])) - *ti); //matches Matlab: hmin=16*eps(t)

    //estimate first step size from derivative at initial condition
    if (*ti == tspan[0])
    {
        //newDt=RCONST(10.0)*hmin;
        //newDt=fmin(sp->dtmax, fabs(tspan[1]-tspan[0]));
        for (int j = 0; j < N_VAR; ++j)
            err[j] = k1[j] / fmax(fabs(xi[j]), threshold);
        
        realtype rh = norm_inf(err, N_VAR) / (RCONST(0.8) * pow(sp->reltol, expon));
        if (newDt * rh > RCONST(1.0))
            newDt = RCONST(1.0) / rh;
            
        newDt=clamp(newDt, hmin, sp->dtmax);
    }

    bool noFailedSteps = true;
    while (true)
    {
        tNew = *ti;
        for (int j = 0; j < N_VAR; j++)
        {
            newxi[j] = xi[j];
            newk1[j] = k1[j];
        }

        newDt = adaptiveOneStep(&tNew, newxi, newk1, pars, newDt, aux, err, wi);
        //returns purified dt: roundoff reduces accuracy of ti+dt, so use the portion of dt that had an effect...

        //Error estimation - elementwise
        for (int j = 0; j < N_VAR; j++)
            err[j] /= fmax( fmax( fabs(xi[j]), fabs(newxi[j]) ) , threshold);
        

        normErr = norm_inf(err, N_VAR); //largest relative error among variables
        // normErr = norm_2(err, N_VAR);
        // normErr = norm_1(err, N_VAR);

        //Error estimation - normcontrol
        //~ realtype nXi=norm_2(xi, N_VAR);
        //~ realtype nNewXi=norm_2(newxi, N_VAR);
        //~ normErr = norm_2(err, N_VAR)/fmax(fmax(nXi,nNewXi),threshold);

        //shrink dt if too much error
        if (normErr > sp->reltol)
        {
            if (newDt < hmin)
            {
                *dt = hmin;
                return -1;
            } //pass error signal back to main loop..

            if (noFailedSteps)
            { //first failure: shrink proportional to error
                noFailedSteps = false;
                newDt *= fmax(ADAPTIVE_STEP_MAX_SHRINK, RCONST(0.9) * pow(sp->reltol / normErr, expon)); //limited shrink rate
            }
            else
            { //repeated failed step: cut stepsize in half
                newDt *= RCONST(0.5) ; 
            }
            
            newDt = clamp(newDt, hmin, sp->dtmax); //limiters
        }
        else
        {
            break;
        }
    }

    //update the solution
    *ti = tNew;
    for (int j = 0; j < N_VAR; j++)
    {
        xi[j] = newxi[j];
        k1[j] = newk1[j];
    }

    //no failure this step => attempt to increase dt for next timestep
    if (noFailedSteps)
    {
        realtype temp = RCONST(0.9) / pow(normErr / sp->reltol, expon); //in case normErr=0
        if (temp < RCONST(5.0))
            newDt *= temp;
        else
            newDt *= RCONST(5.0); //max increase is 5-fold
    }

    newDt = fmin(newDt, fabs(tspan[1] - *ti)); //hit the final time exactly
    newDt = clamp(newDt, hmin, sp->dtmax); //limiters
    
    *dt = newDt; //new step size to attempt on next step

    return 0;
}
#endif

#endif //STEPPERS_H_