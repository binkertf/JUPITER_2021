#include "jupiter.h"

////////////////////////////////////
////////////////////////////////////
////////////////////////////////////
///
/// Isothermal solvers below
///
////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

real ItereRS (pguess, rhoL, rhoR, du, aL, aR)
     real pguess, rhoL, rhoR, du, aL, aR;
{
  real f, fprime;
  f = -aL*log(pguess/rhoL)-aR*(pguess-rhoR)/sqrt(rhoR*pguess)-du;
  fprime = aL/sqrt(pguess)+.5*aR*(1./sqrt(rhoR)+sqrt(rhoR)/pguess);
  fprime = -1./sqrt(pguess)*fprime;
  return -f/fprime + pguess;
}


void GetStar_ITERATIVE (rhoL, rhoR, uL, uR, aL, aR, us, ps)
real rhoL, rhoR, uL, uR, aL, aR;
real *us, *ps;
{
  real deltau, delta, prs, ust;
  real srhoL, srhoR, denom;	/* below: exact 2 shocks */
  rhoL = rhoL*aL*aL;		/* rhoL and rhoR */
  rhoR = rhoR*aR*aR;		/* are now pressures */
  srhoL = sqrt(rhoL);
  srhoR = sqrt(rhoR);
  denom = (aR*srhoL+aL*srhoR)/(srhoL*srhoR);
  deltau = uR-uL;		
  delta = deltau*deltau+4.*(aL*srhoL+aR*srhoR)*denom;
  prs = (-deltau+sqrt(delta))/2./denom;
  prs = prs*prs;
  ust = uR+(prs-rhoR)/sqrt(prs*rhoR)*aR;
  if ((prs < rhoL) && (prs < rhoR)) { /* exact 2 rarefactions */
    prs = pow(rhoL,aL/(aR+aL)) * pow(rhoR,aR/(aR+aL)) * exp(-deltau/(aR+aL));
    ust = uR + aR * log(prs/rhoR);
  }
  if ((prs < rhoL) && (prs > rhoR)) { /* Left rarefaction & right shock */
    prs = ItereRS(prs, rhoL, rhoR, deltau, aL, aR);
    ust = uR+(prs-rhoR)/sqrt(prs*rhoR)*aR;
  }
  if ((prs > rhoL) && (prs < rhoR)) { /* Right rarefaction & left shock */
    prs = ItereRS(prs, rhoR, rhoL, deltau, aR, aL);
    ust = uL+(rhoL-prs)/sqrt(prs*rhoL)*aL;
  }
  *ps = prs;
  *us = ust;
}

void GetStar_TWORAREFACTIONS (rhoL, rhoR, uL, uR, aL, aR, us, ps)
real rhoL, rhoR, uL, uR, aL, aR;
real *us, *ps;
{
  real deltau, prs, ust;
  deltau = uR-uL;		/* below: exact two rarefactions */
  rhoL = rhoL*aL*aL;		/* rhoL and rhoR */
  rhoR = rhoR*aR*aR;		/* are now pressures */
  prs = pow(rhoL,aL/(aR+aL)) * pow(rhoR,aR/(aR+aL)) * exp(-deltau/(aR+aL));
  ust = uR + aR * log(prs/rhoR);
  *ps = prs;
  *us = ust;
}

void GetStar_TWOSHOCKS (rhoL, rhoR, uL, uR, aL, aR, us, ps)
real rhoL, rhoR, uL, uR, aL, aR;
real *us, *ps;
{
  real deltau, delta, prs, ust;
  real srhoL, srhoR, denom;	/* below: exact 2 shocks */
  rhoL = rhoL*aL*aL;		/* rhoL and rhoR */
  rhoR = rhoR*aR*aR;		/* are now pressures */
  srhoL = sqrt(rhoL);
  srhoR = sqrt(rhoR);
  denom = (aR*srhoL+aL*srhoR)/(srhoL*srhoR);
  deltau = uR-uL;		
  delta = deltau*deltau+4.*(aL*srhoL+aR*srhoR)*denom;
  prs = (-deltau+sqrt(delta))/2./denom;
  prs = prs*prs;
  ust = uR+(prs-rhoR)/sqrt(prs*rhoR)*aR;
  *ps = prs;
  *us = ust;
}

////////////////////////////////////
////////////////////////////////////
////////////////////////////////////
///
/// Adiabatic solvers below
///
////////////////////////////////////
////////////////////////////////////
////////////////////////////////////

boolean GetStar_ADI_SS (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR)
     real rhoL, rhoR, uL, uR, aL, aR, gammaL, gammaR;
     real *us, *ps;
{
  real AL, AR, BL, BR;
  real pL=rhoL*aL*aL/gammaL, pR=rhoR*aR*aR/gammaR;
  real df, deltau, a, prs, prs0;
  int count=0;
  deltau = uL - uR;
  prs0 = prs = *ps; // The initial guess (PVRS) is evaluated in the calling function
  AL = (2.0/(gammaL+1.0))/rhoL; 
  AR = (2.0/(gammaR+1.0))/rhoR;
  BL = (gammaL-1.0)/(gammaL+1.0) * pL;
  BR = (gammaR-1.0)/(gammaR+1.0) * pR;
  df = (prs - pL)*sqrt(AL/(prs+BL)) + (prs - pR)*sqrt(AR/(prs+BR));
  while ((fabs(df-deltau)>=1.0e-10) && (count < 20)) {
    count++;
    a = sqrt(AL/(prs+BL))*(1.-0.5*(prs-pL)/(prs+BL)) +		\
      sqrt(AR/(prs+BR))*(1.-0.5*(prs-pR)/(prs+BR));
    prs += -1.0*(df-deltau)/a;
    if (prs < 0.0) prs = 1e-8*prs0;
    //The above corresponds to a too large initial guess (in which
    //case the first pressure iterate is negative). We then take an
    //almost vanishing new initial guess. We could do that always,
    //but generally the PVRS guess (see Toro) needs less
    //iterations.
    df = (prs - pL)*sqrt(AL/(prs+BL)) + (prs - pR)*sqrt(AR/(prs+BR));
  }
  CHECK_CONVERGENCE;
  if ((count >= 20) || (prs < 0.0)) {
    printf ("TSRS convergence failed. Code crashed.");
    INSPECT_REAL (prs);
    INSPECT_REAL (rhoL);
    INSPECT_REAL (rhoR);
    INSPECT_REAL (uL);
    INSPECT_REAL (uR);
    INSPECT_REAL (aL);
    INSPECT_REAL (aR);
    INSPECT_INT (count);
    exit (EXIT_FAILURE);
  }
  *ps = prs;
  *us = 0.5*(uL+uR)+0.5*((prs-pR)*sqrt(AR/(prs+BR))-(prs-pL)*sqrt(AL/(prs+BL))); 
  return NO;
}

boolean GetStar_ADI_SR (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR)
     real rhoL, rhoR, uL, uR, aL, aR, gammaL, gammaR;
     real *us, *ps;
{	/* Left shock, right rarefaction */
  real AL, BL;
  real pL=rhoL*aL*aL/gammaL, pR=rhoR*aR*aR/gammaR; // a is sound speed (Ideal gas)
  real df, deltau, a, prs;
  prs = *ps;  // The initial guess (PVRS) is evaluated in the calling function
  AL = (2.0/(gammaL+1.0))/rhoL;
  BL = (gammaL-1.0)/(gammaL+1.0) * pL;
  deltau = uL - uR;
  df = (prs-pL)*sqrt(AL/(prs+BL)) + 2.0/(gammaR-1.0)*aR*(pow(prs/pR,.5*(gammaR-1.)/gammaR)-1.);
  while (fabs(df-deltau)>=1.e-10){
    a = sqrt(AL/(prs+BL))*(1.-0.5*(prs-pL)/(prs+BL)) + aR/(gammaR*pR)*pow(prs/pR,-.5*(gammaR+1.)/gammaR);
    prs += -(df-deltau)/a;
    if (prs < 0.0) prs = rhoL*aL*aL/gammaL; // Overshoot with initial guess: we adopt the smallest pressure (here, left side) as a new initial guess
    df = (prs-pL)*sqrt(AL/(prs+BL)) + 2.0/(gammaR-1.0)*aR*(pow(prs/pR,.5*(gammaR-1.)/gammaR)-1.);
  }
  CHECK_CONVERGENCE;
  *ps = prs;
  *us = 0.5*( uL + uR )							\
    + 0.5*(2.0/(gammaR-1.0)*aR*(pow(prs/pR,.5*(gammaR-1.)/gammaR)-1.) - (prs-pL)*sqrt(AL/(prs+BL)));
  return NO;
}

boolean GetStar_ADI_RS (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR)
     real rhoL, rhoR, uL, uR, aL, aR, gammaL, gammaR;
     real *us, *ps;
{	/* Right shock, left rarefaction */
  real AR, BR;
  real pL=rhoL*aL*aL/gammaL, pR=rhoR*aR*aR/gammaR;
  real df, deltau, a, prs;
  prs = *ps; // The initial guess (PVRS) is evaluated in the calling function
  AR = (2.0/(gammaR+1.0))/rhoR;
  BR = (gammaR-1.0)/(gammaR+1.0) * pR;
  deltau = uL - uR;
  df = (prs-pR)*sqrt(AR/(prs+BR)) + 2.0/(gammaL-1.0)*aL*(pow(prs/pL,.5/gammaL*(gammaL-1.))-1.);
  while (fabs(df-deltau)>=1.e-10){
    a = sqrt(AR/(prs+BR)) * (1.-0.5*(prs-pR)/(prs+BR)) + aL/(gammaL*pL)*pow(prs/pL,-.5/gammaL*(gammaL+1.));
    prs += -(df-deltau)/a;
    if (prs < 0.0) prs = rhoR*aR*aR/gammaR; // Overshoot with initial guess: we adopt the smallest pressure (here, right side) as a new initial guess
    df = (prs-pR)*sqrt(AR/(prs+BR)) + 2.0/(gammaL-1.0)*aL*(pow(prs/pL,.5/gammaL*(gammaL-1.))-1.);
  }
  CHECK_CONVERGENCE;
  *ps = prs;
  *us = 0.5*( uL + uR )							\
    + 0.5*(-2.0/(gammaL-1.0)*aL*(pow(prs/pL,.5/gammaL*(gammaL-1.))-1.) + (prs-pR)*sqrt(AR/(prs+BR)));
  return NO;
}

boolean GetStar_AdiabaticSolver (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR)
real rhoL, rhoR, uL, uR, aL, aR, gammaL, gammaR;
real *us, *ps;
{
  real egL, egR, prs, ust, rhomean, amean, gammamean;
  real pL,pR;	
  real deltau = uL-uR;

  
  pL = rhoL*aL*aL/gammaL; 
  pR = rhoR*aR*aR/gammaR; 


  rhomean = .5*(rhoL+rhoR);
  amean   = .5*(aL+aR);
  gammamean = .5*(gammaL + gammaR);
  *ps     = .5*(pL+pR)+.5*deltau*rhomean*amean;
  *us     = .5*(uL+uR)+.5*(pL-pR)/(rhomean*amean);

  if ((*ps > pL) && (*ps > pR)) // Two-shock solution
    return GetStar_ADI_SS (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR);

  if ((pL > *ps) && (*ps > pR)) // Left Rarefaction, right shock solution
    return GetStar_ADI_RS (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR);
 
  if ((pL < *ps) && (*ps < pR)) // Right Rarefaction, left shock solution
    return GetStar_ADI_SR (rhoL, rhoR, uL, uR, aL, aR, us, ps, gammaL, gammaR);
 
  /// Staying on the 2-rarefaction (RR) solver below

  prs = ( deltau*0.5 * (gammamean-1.) + aL+aR );
  if (prs < 0.0) // <== in that case vacuum is generated
    return YES;
  egL = (gammaL-1.)/(2.*gammaL);
  egR = (gammaR-1.)/(2.*gammaR);
  prs = prs / ( pow(pL,-egL)*aL + pow(pR,-egR)*aR);
  prs = pow(prs, 2.0/(gammamean-1.0) * gammamean);
  CHECK_CONVERGENCE;
  ust = (uL+uR)*0.5 + aR*(pow(prs/pR,egR)-1.)/(gammaR-1.0) - aL*(pow(prs/pL,egL)-1.) /(gammaL-1.0);
  *ps = prs;
  *us = ust;
  return NO;
}


void GetStar_PRESSURELESS (rhoL, rhoR, uL, uR, aL, aR, us, ps) /* Riemann solver for pressureless fluid */
real rhoL, rhoR, uL, uR, aL, aR;
real *us, *ps;
{
*ps = 0.0; /* pressureless fluid */
*us = (sqrt(rhoL)*uL+sqrt(rhoR)*uR)/(sqrt(rhoL)+sqrt(rhoR)); /* velocity of the delta shock according to Leveque 2004 */
}
