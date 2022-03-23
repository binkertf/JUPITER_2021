#include "jupiter.h"
#include "calc_gamma.h"

#define COOLING !Isothermal


real GetGamma(){
    if(Isothermal) {
        return GAMMA;
    }
    else {
        return GAMMA;
    }
}



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
    double rhs1=0.0, rhs2=0.0, rhs3=0.0;
    double b=0.0,c=0.0, scrh=0.0;
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
    double f[4]={};

    GetSahaHFracs(T, rho, f);
#if HELIUM_IONIZATION == YES
    *mu = 4.0/(  2.*H_MASS_FRAC*(1.0 + f[DEG_y] + 2.0*f[DEG_x]*f[DEG_y])
                 + He_MASS_FRAC*(1.0 + f[DEG_z1]*(1.0 + f[DEG_z2])));
#else
    *mu = 4.0/(2.*H_MASS_FRAC*(1.0 + f[DEG_y] + 2.0*f[DEG_x]*f[DEG_y])
             + He_MASS_FRAC);
#endif
}


/* ***************************************************************** */
double InternalEnergyFunc(double T, double rho)
/*!
 *  Compute the gas internal energy as a function of temperature
 *  and fractions (or density)
 *
 *  \param [in]   T   Gas temperature
 *  \param [in]   rho   Density
 *
 *  \return The gas internal energy (\c rhoe) in code units.
 ******************************************************************* */
{
    double eH=0.0, eHe=0.0, eHpH=0.0, eHplus=0.0, eH2=0.0, rhoe=0.0, eHeplus=0.0, eHeplusplus=0.0;
    double func_zetaR=0.0;
    double f[4]={};

    GetSahaHFracs(T, rho, f);

    func_zetaR = 1.5;   /*   to recover ideal EoS  */
    /*   GetFuncDum(T, &func_zetaR); */
    /* = (1.5 + e(rot) + e(vib)), need to implement */

/* -- Estimate contributions to Egas -- */

    eH   = 1.5*H_MASS_FRAC*(1.0 + f[DEG_x])*f[DEG_y];
    eHe  = 3.0*He_MASS_FRAC/8.0;
    eH2  = 0.5*H_MASS_FRAC*(1.0 - f[DEG_y])*func_zetaR;

/* -- constant terms -- */

#if COOLING == NO
    eHpH   = 4.48*CONST_eV*H_MASS_FRAC*f[DEG_y]/(2.0*BOLTZ*T);
    eHplus = 13.60*CONST_eV*H_MASS_FRAC*f[DEG_x]*f[DEG_y]/(BOLTZ*T);

#if HELIUM_IONIZATION == YES
    eHeplus = 24.59*CONST_eV*He_MASS_FRAC*f[DEG_z1]*(1.0 - f[DEG_z2])/(4.0*CONST_kB*T);
    eHeplusplus = 54.42*CONST_eV*He_MASS_FRAC*f[DEG_z1]*f[DEG_z2]/(4.0*CONST_kB*T);
#else
    eHeplus = eHeplusplus = 0.0;
#endif

#else
    eHpH = eHplus = eHeplus = eHeplusplus = 0.0;
#endif

/* ----------------------------------------------------------------
    Compute rhoe in cgs except for density which is left in
    code units (to avoid extra multiplication and division later)
   ---------------------------------------------------------------- */

    rhoe = (eH2 + eH + eHe + eHpH + eHplus + eHeplus + eHeplusplus)*(BOLTZ*T*rho/XMH);

    return rhoe/(V0*V0);  /* convert to code units */
}


/*
 * need input: v (primitive variables) = (temperature, rho)
 */
/* ********************************************************************* */

double Gamma1(double temperature, double density)
{
    double gmm1=0.0, pressure=0.0;
    double T=0.0;
    double cv=0.0, mu=0.0, chirho = 1.0, chiT = 1.0, rho=0.0, rhoe=0.0;
    double ep=0.0, em=0.0, delta = 1.e-2;
    double Tp=0.0, Tm=0.0, mup=0.0, mum=0.0, rhop=0.0, rhom=0.0;
    double dmu_dT=0.0, dmu_drho=0.0, de_dT=0.0;

/* ---------------------------------------------
    Obtain pressure and fractions.
   --------------------------------------------- */
    T = temperature; /* temperature*/
    rho = density; /* density*/

    GetMu(T, rho, &mu);
    pressure = T*rho/(KELVIN*mu);


/* ---------------------------------------------------
    Compute cV (Specific heat capacity for constant
    volume) using centered derivative.
    cV will be in code units. The corresponding cgs
    units will be the same velocity^2/Kelvin.
   --------------------------------------------------- */

    Tp = T*(1.0 + delta);
    Tm = T*(1.0 - delta);
    em = InternalEnergyFunc(Tm, rho)/rho; /* in code units */
    ep = InternalEnergyFunc(Tp, rho)/rho; /* in code units */

    de_dT = (ep - em)/(2.0*delta*T);
    cv    = de_dT;  /* this is code units. */


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
    gmm1   = pressure/(cv*T*rho)*chiT*chiT  + chirho;

    return gmm1;
}

