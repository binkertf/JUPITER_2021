#include "jupiter.h"

void plm_adiab (beam, dt)
     Beam *beam;
     real dt;
{
  long i, j, k, n;
  real rhoL, rhoR, a, uL, uR, vperpL[2], vperpR[2];
  real *dxe[3], dxp, dxm, sa1, sa2, roa, aor, source_vperp[2];
  /* In the expression below, the first index in the dxrho and dxv
     arrays is for the dimension along which the derivative is to be
     taken, while the second in dxv is for which component of the
     velocity one considers (dxv[i][j] = "dv^j/dx^i") */
  real *vperp[2], *dxrho[3], *dxv[3][3], *rho, *u, *e, *inter;
  real rho0, u0, da1p, da2p, da1m, s_ent, dentp, dentm;
  real e0, energ_transv, eL, eR, da2m, cont_transv, euler_transv;
  n   = beam->length;
  u   = beam->u;
  rho = beam->rho;
  e   = beam->cs;		/* volumic internal energy despite of its name */
  inter= beam->intersurface;
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
    e0   =  e[i];      // Volumic internal energy
    rho0 = rho[i];
    a    = sqrt(GAMMA*(GAMMA-1.0)*e0/rho0); // Sound speed
    aor  = a/rho0;
    roa  = rho0/a;
    if (NDIM>1) vperpR[0] = vperpL[0] = vperp[0][i];
    if (NDIM>2) vperpR[1] = vperpL[1] = vperp[1][i];
    cont_transv  = 0.0;		/* Perpendicular terms for continuity eq. */
    euler_transv = 0.0;		/* Perpendicular terms for Euler eq. */
    energ_transv = 0.0;		/* Perpendicular terms for Energy eq. */
    for (k = 0; k < NDIM-1; k++) {
      cont_transv += vperp[k][i]*dxrho[k+1][i]+rho0*dxv[k+1][k+1][i];
      euler_transv += vperp[k][i]*dxv[k+1][0][i];
      energ_transv += vperp[k][i]*dxe[k+1][i]+e0*GAMMA*dxv[k+1][k+1][i];
    }
    
    rhoL = rhoR = rho0;  //   R : right state on zone left interface
    uL   = uR   = u0;    //   L : left state on zone right interface
    eL   = eR   = e0;    //
    
    dxp  = 2.*(beam->edge[i+1]-beam->center[i]); // Mesh may not be even
    dxm  = 2.*(beam->center[i]-beam->edge[i]);   // in which case dxm != dxp
    
    rhoL += -.5*dt*cont_transv;	
    rhoR += -.5*dt*cont_transv;
    
    uL   += -.5*dt*euler_transv; //   .5*dt : half timestep prediction
    uR   += -.5*dt*euler_transv;
    
    eL   += -.5*dt*energ_transv;
    eR   += -.5*dt*energ_transv;
    
    sa1  = (dxe[0][i]+(GAMMA*e0/a)*dxv[0][0][i]);  // Slopes homogeneous
    sa2  = (dxe[0][i]-(GAMMA*e0/a)*dxv[0][0][i]);  // to the energy
    da1p = .25*(dxp-(u0+a)*dt)*sa1;
    da2p = .25*(dxp-(u0-a)*dt)*sa2;
    da1m = .25*(-dxm-(u0+a)*dt)*sa1;
    da2m = .25*(-dxm-(u0-a)*dt)*sa2;
    
    s_ent= (dxrho[0][i]-rho0/(GAMMA*e0)*dxe[0][i]);
    dentp= .5*(dxp-u0*dt)*s_ent;
    dentm= .5*(-dxm-u0*dt)*s_ent;
    
    if ((u0+a) >= 0.0) {               /* These tests have to be written like that */
      eL   += da1p;                    /* They must not use e.g. if (u < 0.0) ... else ... */
      rhoL += da1p*rho0/(GAMMA*e0);    /* as the result violates right/left symmetry */     
      uL   += da1p*a/(GAMMA*e0);       /* As a consequence a Keplerian disk with equatorial */      
    }                                  /* symmetry would violate this symmetry */
    if ((u0+a) <= 0.0) {
      eR   += da1m;
      rhoR += da1m*rho0/(GAMMA*e0);
      uR   += da1m*a/(GAMMA*e0);
    }
    if ((u0-a) >= 0.0) {
      eL   += da2p;
      rhoL += da2p*rho0/(GAMMA*e0);
      uL   += -da2p*a/(GAMMA*e0);
    }
    if ((u0-a) <= 0.0) {
      eR   += da2m;
      rhoR += da2m*rho0/(GAMMA*e0);
      uR   += -da2m*a/(GAMMA*e0);
    }
    if (u0 >= 0.0) {
      rhoL += dentp;
    }
    if (u0 <= 0.0) {
      rhoR += dentm;
    }
    for (k = 0; k < NDIM-1; k++) {
      source_vperp[k] = -(GAMMA-1.0)/rho0*dxe[k+1][i]+beam->source_perp[k][i];
      for (j = 0; j < NDIM-1; j++)
	source_vperp[k] += -vperp[j][i]*dxv[j+1][k+1][i];
      
      if (u0 >= 0.0) {
        vperpL[k] += .5*dxv[0][k+1][i]*(dxp-u0*dt);
      }
      if (u0 <= 0.0) {
        vperpR[k] += .5*dxv[0][k+1][i]*(-dxm-u0*dt);
      }
      if (beam->MomCorr[k]) {
        vperpL[k] *= beam->MomCorrRatioL[k][i];
        vperpR[k] *= beam->MomCorrRatioR[k][i];
      }
      vperpL[k] += source_vperp[k]*.5*dt;
      vperpR[k] += source_vperp[k]*.5*dt;
      beam->v_perp_L[k][i+1] = vperpL[k];
      beam->v_perp_R[k][i]   = vperpR[k];
    }
    
    if (rhoL <= _Small_Rho) rhoL = _Small_Rho;
    if (rhoR <= _Small_Rho) rhoR = _Small_Rho;
    
    beam->eL[i+1] = eL;
    beam->eR[i] = eR;
    beam->rhoL[i+1] = rhoL;
    beam->rhoR[i] = rhoR;
    beam->uL[i+1] = uL;
    beam->uR[i] = uR;
  }
  AdjustBeamBoundaries (beam);
}
