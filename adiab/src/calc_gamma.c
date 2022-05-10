#include "jupiter.h"


#define NONZERO_INITIALIZE YES /* Fill arrays to nonsense values to catch
                                  uninitialized values later in the code */

/* Switch to choose the molecular hydrogen spin states.
* 0 : Only Para hydrogen ,
* 1 : Equilibrium,
* 2 : Ortho to para ratio of 3:1.
*/
int ORTHO_PARA_MODE = 2;

char     *Array1D (int, int);
#define ARRAY_1D(nx,type)          (type    *)Array1D(nx,sizeof(type))


real GetGamma(){
return GAMMA;
}


/* ********************************************************************* */

char *Array1D (int nx, int dsize)
/*!
 * Allocate memory for a 1-D array of any basic data type starting
 * at 0.
 *
 * \param [in] nx    number of elements to be allocated
 * \param [in] dsize data-type of the array to be allocated
 *
 * \return A pointer of type (char *) to the allocated memory area.
 *********************************************************************** */
{
    char *v;
    v = (char *) malloc (nx*dsize);
    if(!v){
        prs_exit (EXIT_FAILURE);
    }


#if NONZERO_INITIALIZE == YES
    if (dsize==sizeof(double)){
        int i;
        double *q;
        q = (double *) v;
        for (i = nx; i--; ) q[i] =sqrt(-1.0);
    }
#endif
    return v;
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
    rhs3  = 4.0*XMH /rho*rhs2*sqrt(rhs2)*exp(-24.59*CONST_eV/kT);
    b     = 4.0/He_MASS_FRAC*(H_MASS_FRAC + rhs3);
    c     = -rhs3*4.0/He_MASS_FRAC;
    fdeg[DEG_z1] = -2.0*c/(b + sqrt(b*b - 4.0*c));

    rhs3  = XMH /rho*rhs2*sqrt(rhs2)*exp(-54.42*CONST_eV/kT);
    b     = 4.0/He_MASS_FRAC*(H_MASS_FRAC + 0.25*He_MASS_FRAC + rhs3);
    c     = -rhs3*4.0/He_MASS_FRAC;
    fdeg[DEG_z2] = -2.0*c/(b + sqrt(b*b - 4.0*c));
#endif
}
/* ********************************************************************* */



void MakeZetaTables(double *lnT, double *funcdum, int nsteps)
/****************************************************************************
 *  \param [in]  lnT    Array of Logarithmic values of Gas
 *                      temperatures from 0.01 K to 10^12 K.
 *  \param [in] nsteps  Number of equal spacings in T.
 *  \param [out] funcdum The function of zetaR that goes into EH2.
 *
 *  \return This function has no return value
 *********************************************************************** */
{
    int    i,j;
    double Temp0 = 0.01*T_CUT_RHOE;
    double Tmax  = 1.0e12;
    double dy = log(Tmax/Temp0)*(1./nsteps);
    double dT = Temp0*exp(dy);
    double T =0., a=0., b=0., b1=0., scrh=0., inv_T2=0.;
    double zetaP=0., dzetaP=0., zetaO=0., dzetaO=0., zetaR=0., dzetaR=0.;
    double dum1=0., dum2=0., dum3=0.;
    double alpha=0., beta=0., gamma=0.;
    double dzO_zO_m=0., db=0., sum1=0., sum2=0.;


    if (ORTHO_PARA_MODE == 0){
        alpha = 1.0; beta = 0.0; gamma = 0.0;
    }else if(ORTHO_PARA_MODE == 2){
        alpha = 0.25; beta = 0.75; gamma = 0.0;
    }else{
        alpha = 1.0; beta = 0.0; gamma = 1.0;
    }

    b1 = 2.0*THETA_R;
    for(j = 0; j < nsteps; j++){
        T = Temp0*exp(j*dy);
        inv_T2 = 1.0/(T*T);
        zetaO = zetaP = dzetaP = dzetaO = 0.0;
        dzO_zO_m = sum1 = sum2 = 0.0;
        for(i = 0; i <= 10000; i++){
            a = 2*i + 1;
            b = i*(i + 1)*THETA_R;
            if (i%2 == 0){
                scrh    = a*exp(-b/T);
                zetaP  += scrh;
                dzetaP += scrh*b;
            }else{
                db    = b - b1;
                scrh  = a*exp(-db/T);
                sum1 += scrh;
                sum2 += scrh*db;
            }
        }
        dzetaP *= inv_T2;

        zetaO  = exp(-b1/T)*sum1;
        dzetaO = exp(-b1/T)*(b1*sum1 + sum2)*inv_T2;

        dzO_zO_m = sum2/sum1*inv_T2; /* = zeta'(O)/zeta(O) - 2*theta/T^2 */

        lnT[j]  = log(T);

        /* -----------------------------------------
            Compute table
           ----------------------------------------- */

        scrh   = zetaO*exp(2.0*THETA_R/T);

        zetaR  = pow(zetaP,alpha)*pow(scrh,beta) + 3.0*gamma*zetaO;
        dzetaR = (zetaR - 3.0*gamma*zetaO)*(alpha*(dzetaP/zetaP) +
                                            beta*dzO_zO_m)  + 3.0*gamma*dzetaO;
        dum1  = THETA_V/T;
        dum2  = dum1*exp(-dum1)/(1.0 - exp(-dum1));
        dum3  = (T/zetaR)*dzetaR;
        funcdum[j] = 1.5 + dum2 + dum3;
    }
}


/* ********************************************************************* */

void GetFuncDum(double T, double *funcdum_val)
/*!
 *  Interpolate the value of a function of \c zetaR from the table
 *  that is used to estimate the value of EH2.
 *
 * \param [in]   T             Value of temperature in kelvin.
 * \param [out] *funcdum_val   Pointer to the value of function of zetaR
 *
 * \return This function has no return value.
 *
 * \b  Reference:\n
 *     D'Angelo, G. et al, ApJ 778 2013.
 *********************************************************************** */
{
    int    klo, khi, kmid;
    int nsteps = 5000;
    double mu=0., Tmid=0., dT=0.;
    static double *lnT, *funcdum;
    double y=0., dy=0.;
    int indx;

/* -------------------------------------------
    Make table on first call
   ------------------------------------------- */

    if (lnT == NULL){
        lnT    = ARRAY_1D(nsteps, double);
        funcdum = ARRAY_1D(nsteps, double);
        MakeZetaTables(lnT, funcdum, nsteps);
    }
    y = log(T);

/* -------------------------------------------------
    Since the table has regular spacing in log T,
    we divide by the increment to find the nearest
    node in the table.
   ------------------------------------------------- */

    if (y > lnT[nsteps-2]) {
        *funcdum_val = funcdum[nsteps-2];
    } else if (y < lnT[0]) {
        *funcdum_val = funcdum[0];
    } else{
        dy   = lnT[1] - lnT[0];
        indx = floor((y - lnT[0])/dy);

        if (indx >= nsteps || indx < 0){
            printf("! GetFuncDum: indx out of range, indx = %d\n",indx);
            printf("! T = %12.6e\n",T);
            prs_exit (EXIT_FAILURE);
        }
        *funcdum_val = (funcdum[indx]*(lnT[indx+1] - y)
                        + funcdum[indx+1]*(y - lnT[indx]))/dy;
    }
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

    /* func_zetaR = 1.5;   to recover ideal EoS  */
    if(T<10){
        func_zetaR = 1.5;
    }
    /* case of changing gamma */
    else {
        GetFuncDum(T, &func_zetaR);
    }


/* -- Estimate contributions to Egas -- */
    eH   = 1.5*H_MASS_FRAC*(1.0 + f[DEG_x])*f[DEG_y];
    eHe  = 3.0*He_MASS_FRAC/8.0;
    eH2  = 0.5*H_MASS_FRAC*(1.0 - f[DEG_y])*func_zetaR;

/* -- constant terms -- */

    eHpH   = 4.48*CONST_eV*H_MASS_FRAC*f[DEG_y]/(2.0*BOLTZ*T);
    eHplus = 13.60*CONST_eV*H_MASS_FRAC*f[DEG_x]*f[DEG_y]/(BOLTZ*T);

#if HELIUM_IONIZATION == YES
    eHeplus = 24.59*CONST_eV*He_MASS_FRAC*f[DEG_z1]*(1.0 - f[DEG_z2])/(4.0*BOLTZ*T);
    eHeplusplus = 54.42*CONST_eV*He_MASS_FRAC*f[DEG_z1]*f[DEG_z2]/(4.0*BOLTZ*T);
#else
    eHeplus = eHeplusplus = 0.0;
#endif


/* ----------------------------------------------------------------
    Compute rhoe in cgs except for density which is left in
    code units (to avoid extra multiplication and division later)
   ---------------------------------------------------------------- */

    rhoe = (eH2 + eH + eHe + eHpH + eHplus + eHeplus + eHeplusplus)*(BOLTZ*T*rho/XMH);

    return rhoe/(V0*V0);  /* convert to code units */
}



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

    if(HELIUM_IONIZATION){
        /* matches d'Angelo et al. paper */
        ORTHO_PARA_MODE = 2;
    }

    if(density>1e-7){
        /* matches d'Angelo et al. paper */
        ORTHO_PARA_MODE = 2;
    }
    else{
        /* matches B. Vaidya et al. paper */
        ORTHO_PARA_MODE = 1;
    }

    T = temperature; /* temperature*/
    rho = density/RHO0; /* density, needs to be normalized */

    GetMu(T, rho, &mu);
    /* printf("%g %g\n", T, mu); */
    pressure = (T*rho)/(KELVIN*mu);

/* ---------------------------------------------------
    Compute cV (Specific heat capacity for constant
    volume) using centered derivative.
    cV will be in code units. The corresponding cgs
    units will be the same velocity^2/Kelvin.
   --------------------------------------------------- */

    Tp = T*(1.0 + delta);
    Tm = T*(1.0 - delta);
    /* printf("%g %g\n", T, InternalEnergyFunc(T, rho)*XMH/(rho*RHO0)); */
    em = InternalEnergyFunc(Tm, rho)/(rho); /* in code units */
    ep = InternalEnergyFunc(Tp, rho)/(rho); /* in code units */

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
    /* printf("%g %g\n", T, cv);*/

    return gmm1;
}

