#include "jupiter.h"

void plm (beam, dt)
     Beam *beam;
     real dt;
{
  long i, j, k, n, dim;
  real rhoL, rhoR, a, uL, uR, vperpL[2], vperpR[2];
  real dxp, dxm, sa1, sa2, roa, aor, source_vperp[2], rhoc=0.0;
  /* In the expression below, the first index in the dxrho and dxv
     arrays is for the dimension along which the derivative is to be
     taken, while the second in dxv is for which component of the
     velocity one considers (dxv[i][j] = "dv^j/dx^i") */
  real *vperp[2], *dxrho[3], *dxv[3][3], *rho, *u, *cs, *inter;
  real rho0, u0, da1p, da2p, da1m, da2m, cont_transv, euler_transv;

  real radius, radius_in, radius_out, vv0;

  boolean HSE;

  dim = beam->dim[0];
  n   = beam->length;
  u   = beam->u;
  rho = beam->rho;
  cs   = beam->cs;
  inter= beam->intersurface;
  HSE = CorrHydroStat[dim];
  for (k = 0; k < NDIM-1; k++)
    vperp[k] = beam->v_perp[k];
  for (k = 0; k < NDIM; k++) {
    dxrho[k] = beam->slope_rho[k];
    dxv[k][0] = beam->slope_u[k];
    dxv[k][1] = beam->slope_v_perp[0][k];
    dxv[k][2] = beam->slope_v_perp[1][k];
  }
  for (i = 1; i < n-1; i++) {	/* Do not change these limits. Note
				   that in case NGHOST>2, they can
				   probably be changed. However this
				   choice is safe whatever the value
				   of NGHOST */
    u0  = u[i];
    a    =  cs[i];
    rho0 = rho[i];
    aor  = a/rho0;
    roa  = rho0/a;
    if (HSE) {
      rhoc = beam->HS_cent_rho[i];
      rho0 -= rhoc;
    }
    if (NDIM>1) vperpR[0] = vperpL[0] = vperp[0][i];
    if (NDIM>2) vperpR[1] = vperpL[1] = vperp[1][i];
    cont_transv = 0.0;		/* Perpendicular terms for continuity eq. */
    for (k = 0; k < NDIM-1; k++)
      cont_transv += vperp[k][i]*dxrho[k+1][i]+rho0*dxv[k+1][k+1][i];
    euler_transv = 0.0;		/* Perpendicular terms for Euler eq. */
    for (k = 0; k < NDIM-1; k++)
      euler_transv += vperp[k][i]*dxv[k+1][0][i];

    rhoL = rhoR = rho0;
    uL   = uR   = u0;

    dxp  = 2.*(beam->edge[i+1]-beam->center[i]);
    dxm  = 2.*(beam->center[i]-beam->edge[i]);

    if (!HSE) {
      rhoL += -.5*dt*cont_transv;
      rhoR += -.5*dt*cont_transv;
    }				/* cuidadin revoir tout ça */
    /* doit-on sortir les deux lignes suivantes du test ? */
    uL += -.5*dt*euler_transv;
    uR += -.5*dt*euler_transv;

    if (!HSE) {
      rhoL += -.5*dt*(beam->srcrho[i]);
      rhoR += -.5*dt*(beam->srcrho[i]);
      uL   += +.5*dt*(beam->source[i]);
      uR   += +.5*dt*(beam->source[i]);
    }

    sa1  = (dxrho[0][i]+roa*dxv[0][0][i]);
    sa2  = (dxrho[0][i]-roa*dxv[0][0][i]);
    da1p = .25*(dxp-(u0+a)*dt)*sa1;
    da2p = .25*(dxp-(u0-a)*dt)*sa2;
    da1m = .25*(-dxm-(u0+a)*dt)*sa1;
    da2m = .25*(-dxm-(u0-a)*dt)*sa2;

    if ((u0+a) >= 0.0) {    /* These tests have to be written like that */
      rhoL += da1p;       /* They must not use '>' or '<' */
      /*[e.g. if (u < 0.0) ... else ...] */
      uL += da1p*aor;	    /* as the result woud violate */
      /* right/left symmetry */
    }
    if ((u0+a) <= 0.0) {    /* As a consequence a Keplerian disk */
      /* with equatorial symmetry */
      rhoR += da1m;       /* would violate this symmetry */
      uR += da1m*aor;
    }
    if ((u0-a) >= 0.0) {
      rhoL += da2p;
      uL += -da2p*aor;
    }
    if ((u0-a) <= 0.0) {
      rhoR += da2m;
      uR += -da2m*aor;
    }
    for (k = 0; k < NDIM-1; k++) {
      source_vperp[k] = -a*aor*dxrho[k+1][i]+beam->source_perp[k][i];
      for (j = 0; j < NDIM-1; j++)
	source_vperp[k] += -vperp[j][i]*dxv[j+1][k+1][i];
      vperpL[k] += .5*dxv[0][k+1][i]*(dxp-u0*dt);
      vperpR[k] += .5*dxv[0][k+1][i]*(-dxm-u0*dt);

      if ((KEPLERIAN) && (dim == _RAD_) && (beam->dim[k+1] == _AZIM_)) {
	radius_out = beam->edge[i+1];
	radius_in  = beam->edge[i];
	radius = beam->center[i];
	vv0 = radius*OMEGAFRAME;
	beam->v_perp_L[k][i+1] =\
	  (vperp[k][i]+vv0)*pow(radius_out/radius,-.5)-radius_out*OMEGAFRAME;
	beam->v_perp_R[k][i] =\
	  (vperp[k][i]+vv0)*pow(radius_in/radius,-.5)-radius_in*OMEGAFRAME;
      } else {

	//	if (beam->MomCorr[k]) {// cuidadin
	//       vperpL[k] *= beam->MomCorrRatioL[k][i];
	//       vperpR[k] *= beam->MomCorrRatioR[k][i];
	//     }
	//            vperpL[k] += source_vperp[k]*.5*dt;
	//            vperpR[k] += source_vperp[k]*.5*dt;
	//beam->v_perp_L[k][i+1] = vperpL[k];
	//beam->v_perp_R[k][i]   = vperpR[k];
	beam->v_perp_L[k][i+1] = vperp[k][i];
	beam->v_perp_R[k][i]   = vperp[k][i]; /* gfo like... cuidadin */
      }
    }

    if (rhoL+rhoc <= _Small_Rho) rhoL = _Small_Rho-rhoc;
    if (rhoR+rhoc <= _Small_Rho) rhoR = _Small_Rho-rhoc;

    /* The correction below takes into account the gradient of the
       sound speed along the beam, as a source term on the velocity */

    if (!HSE) {
      uL += -.5*dt*(cs[i+1]*cs[i+1]-cs[i]*cs[i])/dxp;
      uR += .5*dt*(cs[i-1]*cs[i-1]-cs[i]*cs[i])/dxm;
    }

    beam->rhoL[i+1] = rhoL;
    beam->rhoR[i] = rhoR;
    beam->uL[i+1] = uL;
    beam->uR[i] = uR;
    if (HSE) {
      beam->rhoL[i+1] += beam->HS_int_rho[i+1];
      beam->rhoR[i] += beam->HS_int_rho[i];
    }
  }
  AdjustBeamBoundaries (beam);
}
