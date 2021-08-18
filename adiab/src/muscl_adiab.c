#include "jupiter.h"

void muscl_adiab (beam, dt)
     Beam *beam;
     real dt;
{
  long i, j, k, n;
  real *vperp[2], *dxrho[3], *dxe[3], *dxv[3][3], uL, uR;
  real *rho, *u, *e, dxp, dxm, rhoL, rhoR, eL, eR, min_e, min_rho;
  (void) dt;
/* In the expression above, the first index in the dxrho and dxv
   arrays is for the dimension along which the derivative is to be
   taken, while the second in dxv is for which component of the
   velocity one considers (dxv[i][j] = "dv^j/dx^i") */

  u   = beam->u_pred;
  rho = beam->rho_pred;
  n   = beam->length;
  e   = beam->e_pred;
  for (k = 0; k < NDIM-1; k++)
    vperp[k] = beam->v_perp_pred[k];
  for (k = 0; k < NDIM; k++) {
    dxrho[k] = beam->slope_rho[k];
    dxe[k] = beam->slope_energy[k];
    dxv[k][0] = beam->slope_u[k];
    if (NDIM > 1)
      dxv[k][1] = beam->slope_v_perp[0][k];
    if (NDIM > 2)
      dxv[k][2] = beam->slope_v_perp[1][k];
  }
  for (k = 0; k < NDIM-1; k++)
    vperp[k] = beam->v_perp_pred[k];
  for (i = 1; i < n-1; i++) {	/* Do not change these limits. Note
				   that in case NGHOST>2, they can
				   probably be changed. However this
				   choice is safe whatever the value
				   of NGHOST */
    
    dxp  = beam->edge[i+1]-beam->center[i];
    dxm  = beam->center[i]-beam->edge[i];

    rhoL = rho[i]+dxp*dxrho[0][i];
    rhoR = rho[i]-dxm*dxrho[0][i];
    
    eL   = e[i]+dxp*dxe[0][i];
    eR   = e[i]-dxm*dxe[0][i];

    uL   = u[i]+dxp*dxv[0][0][i];
    uR   = u[i]-dxm*dxv[0][0][i];

    for (j = 0; j < NDIM-1; j++) {
      beam->v_perp_R[j][i]   = vperp[j][i]-dxm*dxv[0][j+1][i];
      beam->v_perp_L[j][i+1] = vperp[j][i]+dxp*dxv[0][j+1][i];
    }

    min_rho = 0.2 * rho[i];
    min_e   = 0.2 * e[i];

    if ((eR <= min_e) || (eL <= min_e) || (rhoL < min_rho) || (rhoR < min_rho)) {
      rhoL = rhoR = rho[i];
      eL = eR = e[i];
      uL = uR = u[i];
      for (j = 0; j < NDIM-1; j++) {
	beam->v_perp_R[j][i]   = vperp[j][i];
	beam->v_perp_L[j][i+1] = vperp[j][i];
      }
    }

    beam->rhoL[i+1] = rhoL;
    beam->rhoR[i]   = rhoR;
    
    beam->eL[i+1]   = eL;
    beam->eR[i]     = eR;

    beam->uL[i+1]   = uL;
    beam->uR[i]     = uR;

  }
  AdjustBeamBoundaries (beam);
}

