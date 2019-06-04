//TODO: Force all steppers to follow FSAL? allows "free" interpolation in any timestep (important for event detection, trajectory output at specified times). if purifying ti (fixed steppers) must get new dxi in main loop.

//TODO: Adaptive stepper seems to fail to keep stepsize small enough to prevent blowups. Something to do with abstol??

//TODO: stepper wishlist: Multi-step method support,  implicit methods

//TODO: does fma(a,b,c)=a*b+c for improved accuracy?
//TODO: fixed stepsize: if always purifying ti in main loop (t0+step*dt), don't update ti here?
//TODO: test speed increases for "minimizing registers" vs reuse of computations (e.g. h2=*dt*0.5)

#ifndef STEPPERS_H_
#define STEPPERS_H_

#include "clODE_struct_defs.h" //for SolverParams struct definition

//forward declaration of the RHS function
void getRHS(realtype t, realtype x_[], realtype p_[], realtype dx_[], realtype aux_[], realtype w_[]);


// FIXED STEPSIZE EXPLICIT METHODS

/* The fixed steppers use stepcount to purify the T values (eliminates roundoff)
 * This means we duplicate the first step's DX computation for trajectory storage, and we need one extra DX computation for the final step
 */

#ifdef FORWARD_EULER
//Forward Euler
void stepper(realtype *ti, realtype xi[], realtype k1[], realtype pars[], const realtype dt, realtype aux[], realtype wi[])
{
    
//compute k1
	getRHS(*ti, xi, pars, k1, aux, wi);

//update to new ti, xi, k1
	for (int k = 0; k < N_VAR; k++) {
		xi[k]=fma(dt, k1[k], xi[k]); }//more accurate?
		//~ xi[k] += dt * k1[k]; }
		
	//~ *ti += dt;
}
#endif


#ifdef HEUN

void stepper(realtype * ti, realtype xi[], realtype k1[], realtype pars[], const realtype dt, realtype aux[], realtype wi[])
{
//Heun's method (improved Euler) time step. Also returns the derivative k1 at (ti,xi), for features

    realtype th=*ti+dt;
    realtype tmp[N_VAR], k2[N_VAR];   

//compute k1
	getRHS(*ti, xi, pars, k1, aux, wi);
    
//compute k2
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(dt, k1[k], xi[k]); }
        //~ tmp[k]=xi[k]+dt*k1[k]; }
        
	getRHS(th, tmp, pars, k2, aux, wi);
	 
//update to new ti, xi, k1
    for (int k=0;k<N_VAR;k++)
        xi[k]+=dt*RCONST(0.5)*(k1[k] + k2[k]);
	
	//*ti = th;
	
}


#endif




#ifdef RUNGE_KUTTA4

//Explicit Runge-Kutta 4 time step. Also returns the derivative k1 at ti,xi, for features.
void stepper(realtype * ti, realtype xi[], realtype k1[], realtype pars[],const realtype dt, realtype aux[], realtype wi[])
{
    realtype h2=dt*RCONST(0.5), th2=*ti+h2, th=*ti+dt;
    realtype tmp[N_VAR], k2[N_VAR], k3[N_VAR], k4[N_VAR];   

//compute k1
    getRHS(th, xi, pars, k1, aux, wi);

//compute k2
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(h2, k1[k], xi[k]); }
        //~ tmp[k]=xi[k]+h2*k1[k];
    getRHS(th2, tmp, pars, k2, aux, wi);
    
//compute k3
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(h2, k2[k], xi[k]); }
        //~ tmp[k]=xi[k]+h2*k2[k];
    getRHS(th2, tmp, pars, k3, aux, wi);
    
//compute k4
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(dt, k3[k], xi[k]); }
        //~ tmp[k]=xi[k]+ dt*k3[k];
    getRHS(th, tmp, pars, k4, aux, wi);
 
//update to new ti and xi
    for (int k=0;k<N_VAR;k++)
        xi[k]+=dt*(k1[k]+RCONST(2.0)*k2[k]+RCONST(2.0)*k3[k]+k4[k])/RCONST(6.0);
    
    //*ti = th;
    
}

#endif



// TODO: FIXED STEP MULTI STEP METHODS
// need to fill the first few steps with a single stepper, so if MULTI_STEP is defined also include Euler or RK4 or something
// #define MULTI_STEP 



// ADAPTIVE STEPSIZE EXPLICIT METHODS

#ifdef HEUN_EULER
#define ADAPTIVE_STEPSIZE
#define ORDERPLUSONE	RCONST(2.0)
#define ADAPTIVE_STEP_MAX_SHRINK  RCONST(0.5)

void adaptiveOneStep(realtype *ti, realtype xi[], realtype k1[], realtype pars[], const realtype dt, realtype aux[], realtype err[], realtype wi[])
{
    realtype th=*ti+dt;
    realtype tmp[N_VAR], k2[N_VAR];   

//compute k1
	getRHS(*ti, xi, pars, k1, aux, wi);
    
//compute k2
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(dt, k1[k], xi[k]); }
        //~ tmp[k]=xi[k]+dt*k1[k]; }
        
	getRHS(th, tmp, pars, k2, aux, wi);
	 
//update to new xi
    for (int k=0;k<N_VAR;k++)
        xi[k]+=dt*RCONST(0.5)*(k1[k] + k2[k]);
        
        
//error estimate in each variable
    for (int k=0;k<N_VAR;k++)
		err[k]=xi[k]-tmp[k];
		
	*ti=th;
}

#endif


#ifdef BOGACKI_SHAMPINE23
#define ADAPTIVE_STEPSIZE
#define FSAL_STEP_PROPERTY
#define ORDERPLUSONE	RCONST(3.0)
#define ADAPTIVE_STEP_MAX_SHRINK RCONST(0.5)

#define B1 RCONST(0.222222222222222)
#define B2 RCONST(0.333333333333333)
#define B3 RCONST(0.444444444444444)

/*
#define _B1 RCONST(0.291666666666666)
#define _B2 RCONST(0.250000000000000)
#define _B3 RCONST(0.333333333333333)
#define _B4 RCONST(0.125000000000000)
*/

// _E1=B1-_B1 directly in error estimate dt*(E dot k) = xnew_B - xnewB
#define _E1 RCONST(-0.069444444444444)
#define _E2 RCONST( 0.083333333333333)
#define _E3 RCONST( 0.111111111111111)
#define _E4 RCONST(-0.125000000000000)

//TODO: doensn't match Matlab exactly...

//this is ode23 in matlab
void adaptiveOneStep(realtype *ti, realtype xi[], realtype k1[], realtype pars[], const realtype dt, realtype aux[], realtype err[], realtype wi[])
{
    realtype h2=dt*RCONST(0.5), h3=dt*RCONST(0.75), th=*ti+dt;
    realtype tmp[N_VAR], k2[N_VAR], k3[N_VAR], k4[N_VAR];   

//expects k1 to be precomputed (FSAL)

//compute k2
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(h2, k1[k], xi[k]); }
        //~ tmp[k]=xi[k]+h2*k1[k];
    for (int k=0;k<N_VAR;k++)
        tmp[k]=xi[k]+h2*k1[k];	
    getRHS(*ti+h2, tmp, pars, k2, aux, wi);
    
//compute k3
    for (int k=0;k<N_VAR;k++) {
		tmp[k]=fma(h3, k2[k], xi[k]); }
        //~ tmp[k]=xi[k]+h3*k2[k];
    getRHS(*ti+h3, tmp, pars, k3, aux, wi);

//update xi
	for (int k=0;k<N_VAR;k++)
	{
		//~ tmp[k]=xi[k];
		xi[k]+=dt*(B1*k1[k] +B2*k2[k] +B3*k3[k]); //third order
	}
        
//compute k4
    getRHS(th, xi, pars, k4, aux, wi);
        
//update error estimate
    for (int k=0;k<N_VAR;k++)
    {
        //~ tmp[k]+=dt*(_B1*k1[k] +_B2*k2[k] +_B3*k3[k] +_B4*k4[k]); //second order
		//~ err[k]=xi[k]-tmp[k];
		err[k]=dt*(_E1*k1[k] + _E2*k2[k] + _E3*k3[k] +_E4*k4[k]);
		k1[k]=k4[k]; //first same as last (FSAL) property
	}
		
	*ti=th;
}
#endif


//~ #ifdef CASH_KARP
//~ #define ADAPTIVE_STEPSIZE
//~ #define ORDERPLUSONE	5
//~ 
//~ #endif


#ifdef DORPRI5
#define ADAPTIVE_STEPSIZE
#define FSAL_STEP_PROPERTY
#define ORDERPLUSONE	RCONST(5.0)
#define ADAPTIVE_STEP_MAX_SHRINK RCONST(0.1)

#define C2 RCONST(0.200000000000000)
#define C3 RCONST(0.300000000000000)
#define C4 RCONST(0.800000000000000)
#define C5 RCONST(0.888888888888889)
#define A21 RCONST(0.200000000000000)
#define A31 RCONST(0.075000000000000)
#define A32 RCONST(0.225000000000000)
#define A41 RCONST(0.977777777777778)
#define A42 RCONST(-3.733333333333333)
#define A43 RCONST(3.555555555555555)
#define A51 RCONST(2.952598689224204)
#define A52 RCONST(-11.595793324188385)
#define A53 RCONST(9.822892851699436)
#define A54 RCONST(-0.290809327846365)
#define A61 RCONST(2.846275252525253)
#define A62 RCONST(-10.757575757575758)
#define A63 RCONST(8.906422717743473)
#define A64 RCONST(0.278409090909091)
#define A65 RCONST(-0.273531303602058)
#define B1 RCONST(0.091145833333333)
#define B3 RCONST(0.449236298292902)
#define B4 RCONST(0.651041666666667)
#define B5 RCONST(-0.322376179245283)
#define B6 RCONST(0.130952380952381)

/*
#define _B1 RCONST(0.089913194444444)
#define _B3 RCONST(0.453489068583408)
#define _B4 RCONST(0.614062500000000)
#define _B5 RCONST(-0.271512382075472)
#define _B6 RCONST(0.089047619047619)
#define _B7 RCONST(0.025000000000000)
*/

//Matlab uses _E1=B1-_B1 directly in error estimate dt*(f dot E)
#define _E1 RCONST(0.001232638888889)
#define _E3 RCONST(-0.004252770290506)
#define _E4 RCONST(0.036979166666667)
#define _E5 RCONST(-0.050863797169811)
#define _E6 RCONST(0.041904761904762)
#define _E7 RCONST(-0.025000000000000)


//NOTE: alternative error coefficients exist, use matlab's?

//this is ode45 in matlab
void adaptiveOneStep(realtype *ti, realtype xi[], realtype k1[], realtype pars[], const realtype dt, realtype aux[], realtype err[], realtype wi[])
{
    realtype th=*ti+dt;
    realtype tmp[N_VAR], k2[N_VAR], k3[N_VAR], k4[N_VAR], k5[N_VAR], k6[N_VAR], k7[N_VAR];   

//expects k1 to be precomputed (FSAL)

//compute k2
    for (int k=0;k<N_VAR;k++)
        tmp[k]=xi[k]+A21*dt*k1[k];
    getRHS(*ti+C2*dt, tmp, pars, k2, aux, wi);
    
//compute k3
    for (int k=0;k<N_VAR;k++)
        tmp[k]=xi[k]+dt*(A31*k1[k]+A32*k2[k]);
    getRHS(*ti+C3*dt, tmp, pars, k3, aux, wi);

//compute k4
    for (int k=0;k<N_VAR;k++)
        tmp[k]=xi[k]+dt*(A41*k1[k]+A42*k2[k]+A43*k3[k]);
    getRHS(*ti+C4*dt, tmp, pars, k4, aux, wi);
    
//compute k5
    for (int k=0;k<N_VAR;k++)
        tmp[k]=xi[k]+dt*(A51*k1[k]+A52*k2[k]+A53*k3[k]+A54*k4[k]);
    getRHS(*ti+C5*dt, tmp, pars, k5, aux, wi);
    
//compute k6
    for (int k=0;k<N_VAR;k++)
        tmp[k]=xi[k]+dt*(A61*k1[k]+A62*k2[k]+A63*k3[k]+A64*k4[k]+A65*k5[k]);
    getRHS(th, tmp, pars, k6, aux, wi);
    
//update xi
	for (int k=0;k<N_VAR;k++)
	{
		//~ tmp[k]=xi[k];
		xi[k]+=dt*(B1*k1[k] +B3*k3[k] +B4*k4[k] +B5*k5[k] +B6*k6[k]); //fifth order
	}
        
//compute k7
    getRHS(th, xi, pars, k7, aux, wi);
        
//update error estimate
    for (int k=0;k<N_VAR;k++)
    {
        //~ tmp[k]+=dt*(_B1*k1[k] +_B3*k3[k] +_B4*k4[k] +_B5*k5[k] +_B6*k6[k] +_B7*k7[k]); //fourth order
		//~ err[k]=xi[k]-tmp[k];
        err[k]=dt*(_E1*k1[k] +_E3*k3[k] +_E4*k4[k] +_E5*k5[k] +_E6*k6[k] +_E7*k7[k]); //fourth order
		k1[k]=k7[k]; //first same as last (FSAL) property
	}
		
	*ti=th;
}

#endif


//TODO: modify to work in my context
// low-register use adaptive stepsize solvers, from rktides
//~ #ifdef RK435_KENNEDY
//~ #define ADAPTIVE_STEPSIZE
//~ #endif

//~ #ifdef RK549_KENNEDY
//~ #define ADAPTIVE_STEPSIZE
//~ #endif



#ifdef ADAPTIVE_STEPSIZE

//Wrapper to handle step-size adaptation.  note: wi should be zeros
int stepper(realtype *ti, realtype xi[], realtype k1[], realtype pars[], __constant struct SolverParams* sp, realtype *dt, __constant realtype tspan[], realtype aux[], realtype wi[])
{
	realtype err[N_VAR];
	realtype normErr, newDt=*dt;
	
	realtype threshold=sp->abstol/sp->reltol;
	realtype expon=RCONST(1.0)/ORDERPLUSONE;
	realtype hmin=RCONST(16.0)*fabs(fabs(nextafter(*ti,tspan[1]))-*ti); //matches Matlab: hmin=16*eps(t)
	
	//estimate first step size from derivative at initial condition
	if (*ti == tspan[0]) {
		newDt = fmin(sp->dtmax, fabs(tspan[1]-tspan[0]));
		for (int j=0; j<N_VAR; ++j) { 
			err[j]=k1[j] / fmax(fabs(xi[j]),threshold);
		}
		realtype rh = norm_inf(err, N_VAR) / ( RCONST(0.8)*pow(sp->reltol, expon) );
		if (newDt*rh > RCONST(1.0)) {
			newDt = RCONST(1.0) / rh;
		}
		newDt = fmax(newDt, hmin);
	}
	
	//limit newDt into [hmin, sp->dtmax]
	newDt=clamp(newDt, hmin, sp->dtmax);  // equivalent to: newDt=fmin(fmax(newDt, hmin), sp->dtmax);
	
	//hit the final time exactly
	newDt=fmin(newDt, fabs(tspan[1]-*ti));
	
	realtype tNew, newxi[N_VAR], newk1[N_VAR];
	
	bool nofailed=true;
	while (true)
	{
		tNew=*ti;
		for (int j=0; j<N_VAR; j++)
		{
			newxi[j] = xi[j];
			newk1[j] = k1[j];
		}
		
		adaptiveOneStep(&tNew, newxi, newk1, pars, newDt, aux, err, wi);
	
		//purify dt: roundoff reduces accuracy of ti+dt, so use the portion of dt that had an effect...
		newDt=tNew-*ti;
	
		//Error estimation - elementwise
		for (int j=0; j<N_VAR; j++) {
			err[j]/=fmax( fmax( fabs(xi[j]), fabs(newxi[j]) ), threshold); }
			
		normErr=norm_inf(err, N_VAR); //largest relative error among variables
		// normErr=norm_2(err, N_VAR);
		// normErr=norm_1(err, N_VAR);
		
		//Error estimation - normcontrol
		//~ realtype nXi=norm_2(xi, N_VAR);
		//~ realtype nNewXi=norm_2(newxi, N_VAR);
		//~ normErr=norm_2(err, N_VAR)/fmax(fmax(nXi,nNewXi),threshold);
		
		//shrink dt if too much error
		if (normErr >= sp->reltol)
		{
			if (nofailed) { //first failure: shrink proportional to error
				nofailed=false;
				newDt*=fmax(ADAPTIVE_STEP_MAX_SHRINK, (RCONST(0.8)*pow(sp->reltol/normErr, expon)) );
			}
			else {
				newDt*=RCONST(0.5);
			}
			newDt=fmax(newDt, hmin);
			
			if (newDt<=hmin) { return -1; } //pass some error signal back to main loop
		}
		else 
		{ break; }
	}
	
	//no failure this step => attempt to increase dt for next timestep
	if (nofailed) {
		realtype temp = RCONST(1.25) * pow(normErr/sp->reltol, expon); //in case normErr=0
		if (temp>RCONST(0.2)) {
			newDt/=temp;
		}
		else {
			newDt*=RCONST(5.0);
		}
	}
	
	//update the solution
	*dt = newDt;
	*ti = tNew;
	for (int j=0; j<N_VAR; j++)
	{
		xi[j] = newxi[j];
		k1[j] = newk1[j];
	}
	
	return 0;
}
#endif


#endif //STEPPERS_H_

/*
//Wrapper to handle step-size adaptation
void stepper(realtype *ti, realtype xi[], realtype k1[], realtype pars[], __constant struct SolverParams* sp, realtype *dt, __constant realtype tspan[], realtype aux[], realtype wi[])
{
	realtype err[N_VAR], newxi[N_VAR], newk1[N_VAR];
	realtype tOld, tNew, dtScale, normErr;
	
	realtype tempWi[N_WIENER];
	realtype threshold=sp->abstol/sp->reltol;
	
	tNew=*ti;
	
	//estimate first step from derivative at initial point
	if (*ti== tspan[0]) {
		*dt = min(sp->dtmax, fabs(tspan[1]-tspan[0]));
		realtype tmp[N_VAR];
		for (int j=0; j<N_VAR; ++j) { 
			tmp[j]=k1[j] / MAX(fabs(xi[j]),threshold);
		}
		realtype rh = norm_inf(tmp, N_VAR) / ( 0.8*pow(sp->reltol, RCONST(1.0)/ORDERPLUSONE) );
		if (*dt*rh > RCONST(1.0)) {
			*dt = RCONST(1.0) / rh;
		}
		*dt = max(*dt, sp->dtmin);
	}
	
	for (int j=0; j<N_VAR; ++j)
	{
		newxi[j] = xi[j];
		newk1[j] = k1[j];
	}
				
	bool stepSuccessful=false;
	
	while  (!stepSuccessful)
	{
		
		for (int j=0; j<N_WIENER; ++j)
			tempWi[j]=wi[j]/sqrt(*dt);
		
		adaptiveOneStep(&tNew, newxi, newk1, pars, *dt, aux, err, tempWi);
		
		for (int j=0; j<N_VAR; j++)
			err[j]/=MAX( MAX(fabs(xi[j]),fabs(newxi[j])), threshold);
			
		//infinity norm of error
		//~ normErr=*dt*norm_1(err, N_VAR);
		//~ normErr=*dt*norm_2(err, N_VAR);
		normErr=*dt*norm_inf(err, N_VAR);
		
		//check error
		if (normErr <= sp->reltol)
		{
			
			*ti = tNew;
			
			for (int j=0; j<N_VAR; j++)
			{
				xi[j] = newxi[j];
				k1[j] = newk1[j];
			}
				
			stepSuccessful=true;
		}

		// update step size
		dtScale = 0.8 * pow (sp->reltol/normErr, RCONST(1.0)/ORDERPLUSONE); // multiplier for dt of the next step
		dtScale = MAX (dtScale, 0.1);					  // minimum multiplier = 1/10
		dtScale = MIN (dtScale, 2.0);					  // maximum multiplier = 2
		*dt *= dtScale;				  
		
		*dt = MIN(*dt, tspan[1]-tNew);  	  // Hit final time exactly
		*dt = MAX(*dt, sp->dtmin);
		*dt = MIN(*dt, sp->dtmax);
	}
}
*/
