#include "jupiter.h"

void gfo (beam, dt)
     Beam *beam;
     real dt;
{
  long k, length;
  real *inter;
  boolean chs;
  (void) dt;
  length = beam->length;
  inter = beam->intersurface;
  chs = CorrHydroStat[beam->dim[0]];
  for (k = 1; k < length; k++) {
    beam->rhoL[k] = beam->rho[k-1];
    if (chs)
      beam->rhoL[k] += beam->HS_int_rho[k]-beam->HS_cent_rho[k-1];
    beam->rhoR[k] = beam->rho[k];
    if (chs)
      beam->rhoR[k] += beam->HS_int_rho[k]-beam->HS_cent_rho[k];
    beam->uL[k] = beam->u[k-1];
    beam->uR[k] = beam->u[k];
    beam->v_perp_L[0][k] = beam->v_perp[0][k-1];
    beam->v_perp_R[0][k] = beam->v_perp[0][k];
    beam->v_perp_L[1][k] = beam->v_perp[1][k-1];
    beam->v_perp_R[1][k] = beam->v_perp[1][k];
  }
  AdjustBeamBoundaries (beam);
}
