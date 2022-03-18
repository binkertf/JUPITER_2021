#include "jupiter.h"
#include "calc_gamma.h"

real GetGamma(){
    if(Isothermal) {
        return GAMMA;
    }
    else {
        return GAMMA;
    }
}



/*
 * need input: v (primitive variables) = (rho, pressure, velocity, temperature?)
 * need functions: GetPV_Temperature, GetMu, MeanMolecularWeight, InternalEnergy, GetHFracs
 * !!! code units
 */



/* ********************************************************************* */
void GetSahaHFracs(double T, double rho, double *fdeg)
/*!
 * Compute degree of ionization and dissociation using Saha Equations.
 * The quadratic equation ay^2 + by + c = 0  is solved using
 * a = 1, c = -b (for hydrogen). A similar expression is used for Helium.
 *
 * \param [in]   T     Gas temperature in Kelvin.
 * \param [in]   rho   Gas density (in code units)
 * \param [out]  fdeg  array of ionization/dissociation degrees.
 *********************************************************************** */
{
    double rhs1, rhs2, rhs3;
    double b,c, scrh;
    double kT   = BOLTZ*T;
    double hbar = CONST_h/(2.0*PI);

    rho *= RHO0;

    rhs1 = XMH/(2.0*H_MASS_FRAC*rho);
    rhs2 = XMH*kT/(4.0*PI*hbar*hbar);
    rhs3 = exp(-4.48*CONST_eV/kT);
    b    = rhs1*rhs2*sqrt(rhs2)*rhs3;

    fdeg[DEG_y] = 2.0/(1.0 + sqrt(1.0 + 4.0/b)); /* solution of quadratic equation */

    rhs1 = XMH/(H_MASS_FRAC*rho);
    rhs2 = CONST_me*kT/(2.0*PI*hbar*hbar);
    rhs3 = exp(-13.60*CONST_eV/kT);
    b    = rhs1*rhs2*sqrt(rhs2)*rhs3;

    fdeg[DEG_x] = 2.0/(1.0 + sqrt(1.0 + 4.0/b)); /* solution of quadratic equation */

#if HELIUM_IONIZATION == YES
    rhs3  = 4.0*CONST_mH/rho*rhs2*sqrt(rhs2)*exp(-24.59*CONST_eV/kT);
    b     = 4.0/He_MASS_FRAC*(H_MASS_FRAC + rhs3);
    c     = -rhs3*4.0/He_MASS_FRAC;
    fdeg[DEG_z1] = -2.0*c/(b + sqrt(b*b - 4.0*c));

    rhs3  = CONST_mH/rho*rhs2*sqrt(rhs2)*exp(-54.42*CONST_eV/kT);
    b     = 4.0/He_MASS_FRAC*(H_MASS_FRAC + 0.25*He_MASS_FRAC + rhs3);
    c     = -rhs3*4.0/He_MASS_FRAC;
    fdeg[DEG_z2] = -2.0*c/(b + sqrt(b*b - 4.0*c));
#endif
}



/* ********************************************************************* */
void GetMu(double T, double rho, double *mu)
/*!
 *  Calculate the mean molecular weight for the case in which
 *  hydrogen fractions are estimated using Saha Equations.
 *
 *  \param [in]   T    Gas temperature in Kelvin.
 *  \param [in]   rho  Gas density (code units)
 *  \param [out]  mu   Mean molecular weight
 *********************************************************************** */
{
    double f[4];

    GetSahaHFracs(T, rho, f);
#if HELIUM_IONIZATION == YES
    *mu = 4.0/(  2.*H_MASS_FRAC*(1.0 + f[DEG_y] + 2.0*f[DEG_x]*f[DEG_y])
                 + He_MASS_FRAC*(1.0 + f[DEG_z1]*(1.0 + f[DEG_z2])));
#else
    *mu = 4.0/(2.*H_MASS_FRAC*(1.0 + f[DEG_y] + 2.0*f[DEG_x]*f[DEG_y])
             + He_MASS_FRAC);
#endif
}



/* ********************************************************************* */

double Gamma1(double *v)
{
    double gmm1, T;
    double cv, mu, chirho = 1.0, chiT = 1.0, rho, rhoe;
    double ep, em, delta = 1.e-2;
    double Tp, Tm, mup, mum, rhop, rhom;
    double dmu_dT, dmu_drho, de_dT;
/* ---------------------------------------------
    Obtain temperature and fractions.
   --------------------------------------------- */
    GetMu(T, v[RHO], &mu);
    /*
     * want to calculate temperature, either:
     *     GetPV_Temperature(v, &T); (lookup table)
     *     T   = v[PRS]/v[RHO]*KELVIN*mu; (direct)
     */


/* ---------------------------------------------------
    Compute cV (Specific heat capacity for constant
    volume) using centered derivative.
    cV will be in code units. The corresponding cgs
    units will be the same velocity^2/Kelvin.
   --------------------------------------------------- */

    Tp = T*(1.0 + delta);
    Tm = T*(1.0 - delta);
    em = InternalEnergy(v, Tm)/v[RHO]; /* in code units */
    ep = InternalEnergy(v, Tp)/v[RHO]; /* in code units */

    de_dT = (ep - em)/(2.0*delta*T);
    cv    = de_dT;  /* this is code units. */

    rho = v[RHO];

    GetMu(Tp, rho, &mup);
    GetMu(Tm, rho, &mum);
    dmu_dT = (mup - mum)/(2.0*delta*T);
    chiT  = 1.0 - T/mu*dmu_dT;

    rhop = rho*(1.0 + delta);
    rhom = rho*(1.0 - delta);
    GetMu(T, rhop, &mup);
    GetMu(T, rhom, &mum);
    dmu_drho = (mup - mum)/(2.0*delta*rho);
    chirho   = 1.0 - rho/mu*dmu_drho;

/* --------------------------------------------
    Compute first adiabatic index
   -------------------------------------------- */
    gmm1   = v[PRS]/(cv*T*v[RHO])*chiT*chiT  + chirho;

    return gmm1;
}
