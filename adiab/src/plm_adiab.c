#include "jupiter.h"

void plm_adiab (beam, dt)
     Beam *beam;
     real dt;
{
  long i, j, k, n, dim, lev;
  real rhoL, rhoR, p0, gamma, a, uL, uR, vperpL[2], vperpR[2];
  real *dxe[3], dxp, dxm, sa1, sa2, source_vperp[2];
  /* In the expression below, the first index in the dxrho and dxv
     arrays is for the dimension along which the derivative is to be
     taken, while the second in dxv is for which component of the
     velocity one considers (dxv[i][j] = "dv^j/dx^i") */
  real *vperp[2], *dxrho[3], *dxv[3][3], *rho, *u, *e, *inter;
  real rho0, u0, da1p, da2p, da1m, s_ent, dentp, dentm;
  real e0, energ_transv, eL, eR, da2m, cont_transv, euler_transv;
  real min_rho, min_e;

  lev = beam->desc->level;
  dim = beam->dim[0];
  n   = beam->length;
  u   = beam->u;
  rho = beam->rho;
  e   = beam->cs;		/* volumic internal energy */
  gamma = GAMMA;
  inter= beam->intersurface;
  if (CorrHydroStat[dim]) 
    prs_error ("PLM adiabatic method incompatible with hydrostatic equilibrium enforcement\n");
  for (k = 0; k < NDIM-1; k++)
    vperp[k] = beam->v_perp[k];
  for (k = 0; k < NDIM; k++) {
    dxrho[k]  = beam->slope_rho[k];
    dxe[k]    = beam->slope_energy[k];
    dxv[k][0] = beam->slope_u[k];
    dxv[k][1] = beam->slope_v_perp[0][k];
    dxv[k][2] = beam->slope_v_perp[1][k];
  }

  for (i = 1; i < n-1; i++) {	/* Do not change these limits. */
    u0   = u[i];
    e0   = e[i];      // Volumic internal energy
    rho0 = rho[i];
    eL   = eR   = e0;    //
    rhoL = rhoR = rho0;  //   R : right state on zone left interface

    if(!Isothermal && Stellar && !CONST_GAMMA){
      p0 = secant_method(e0, rho0, 0.00001); // pressure in code units
      gamma = Gamma1(p0 * MU / rho0 * TEMP0, rho0);
    }
    else{
      gamma = GAMMA;
      p0 = e0*(gamma-1.); // pressure
    }

    a    = sqrt(gamma*(gamma-1.0)*e0/rho0); // Adiabatic sound speed
    if (NDIM>1) vperpR[0] = vperpL[0] = vperp[0][i];
    if (NDIM>2) vperpR[1] = vperpL[1] = vperp[1][i];
    cont_transv  = 0.0;		/* Perpendicular terms for continuity eq. */
    euler_transv = 0.0;		/* Perpendicular terms for Euler eq. */
    energ_transv = 0.0;		/* Perpendicular terms for energy eq. */
    for (k = 0; k < NDIM-1; k++) {
      cont_transv += vperp[k][i]*dxrho[k+1][i] + rho0*dxv[k+1][k+1][i];
      euler_transv += vperp[k][i]*dxv[k+1][0][i];
      energ_transv += vperp[k][i]*dxe[k+1][i] + e0*gamma*dxv[k+1][k+1][i];
    }

    uL   = uR   = u0;    //   L : left state on zone right interface

    dxp  = 2.*(beam->edge[i+1]-beam->center[i]); // Mesh may not be even
    dxm  = 2.*(beam->center[i]-beam->edge[i]);   // in which case dxm != dxp
    
    rhoL += -.5*dt*cont_transv;	
    rhoR += -.5*dt*cont_transv;
    eL   += -.5*dt*energ_transv;
    eR   += -.5*dt*energ_transv;

    uL   += -.5*dt*euler_transv; //   .5*dt : prediction at half time step
    uR   += -.5*dt*euler_transv;

    uL   += .5*dt*beam->source[i];
    uR   += .5*dt*beam->source[i];

    sa1  = dxe[0][i]+(gamma*e0/a)*dxv[0][0][i];  // Slopes homogeneous
    sa2  = dxe[0][i]-(gamma*e0/a)*dxv[0][0][i];  // to that of the energy
    da1p = .25*(dxp-(u0+a)*dt)*sa1;
    da2p = .25*(dxp-(u0-a)*dt)*sa2;
    da1m = .25*(-dxm-(u0+a)*dt)*sa1;
    da2m = .25*(-dxm-(u0-a)*dt)*sa2;

    s_ent= dxrho[0][i]-rho0/(gamma*e0)*dxe[0][i]; // All these quantities are
    dentp= .5*(dxp-u0*dt)*s_ent;                    // homogeneous to density derivatives
    dentm= .5*(-dxm-u0*dt)*s_ent;                   // despite of their names

    if ((u0+a) >= 0.0) {               /* These tests have to be written like that */
      eL   += da1p;                    /* In order to enforce left/right symmetry */
      rhoL += da1p*rho0/(gamma*e0);
      uL   += da1p*a/(gamma*e0);
    }
    if ((u0+a) <= 0.0) {
      eR   += da1m;
      rhoR += da1m*rho0/(gamma*e0);
      uR   += da1m*a/(gamma*e0);
    }
    if ((u0-a) >= 0.0) {
      eL   += da2p;
      rhoL += da2p*rho0/(gamma*e0);
      uL   += -da2p*a/(gamma*e0);
    }
    if ((u0-a) <= 0.0) {
      eR   += da2m;
      rhoR += da2m*rho0/(gamma*e0);
      uR   += -da2m*a/(gamma*e0);
    }
    if (u0 >= 0.0) {
      rhoL += dentp;
    }
    if (u0 <= 0.0) {
      rhoR += dentm;
    }
    for (k = 0; k < NDIM-1; k++) {
      source_vperp[k] = beam->source_perp[k][i];
      for (j = 0; j < NDIM-1; j++)
	source_vperp[k] += -vperp[j][i]*dxv[j+1][k+1][i];

      if (u0 >= 0.0) {
	vperpL[k] += .5*dxv[0][k+1][i]*(dxp-u0*dt);
      }
      if (u0 <= 0.0) {
	vperpR[k] += .5*dxv[0][k+1][i]*(-dxm-u0*dt);
      }
      /* Below: adding the pressure source term on the transversal
	 velocities, as this term is no longer calculated in
	 fillsource_predict.c (if method is PLM) */
      if (beam->dim[k+1] != _COLAT_) {
	vperpL[k] += (source_vperp[k]-(gamma-1.0)*dxe[k+1][i]/rho0)*.5*dt;
	vperpR[k] += (source_vperp[k]-(gamma-1.0)*dxe[k+1][i]/rho0)*.5*dt;
      }
      beam->v_perp_L[k][i+1] = vperpL[k];
      beam->v_perp_R[k][i]   = vperpR[k];
    }

    // Checking for no undershoot in steep parts of the flow.

    min_rho = 0.1 * rho[i];
    min_e = 0.1*e[i];

    if (rhoL <= min_rho) {
      rhoL = min_rho;
      eL = a*a*rhoL/(gamma*(gamma-1.0));
    }

    if (eL <= 0.0) {
      eL = min_e;
      rhoL = eL*gamma*(gamma-1.0)/(a*a);
    }

    if (rhoR <= min_rho) {
      rhoR = min_rho;
      eR = a*a*rhoR/(gamma*(gamma-1.0));
    }

    if (eR <= 0.0) {
      eR = min_e;
      rhoR = eR*gamma*(gamma-1.0)/(a*a);
    }

    beam->eL[i+1] = eL;
    beam->eR[i] = eR;
    beam->rhoL[i+1] = rhoL;
    beam->rhoR[i] = rhoR;
    beam->uL[i+1] = uL;
    beam->uR[i] = uR;

  }
  AdjustBeamBoundaries (beam);
}
