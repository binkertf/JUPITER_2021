#include "jupiter.h"

void gfo_adiab (beam, dt)
     Beam *beam;
     real dt;
{
  long k, length;
  real *inter;
  length = beam->length;
  inter = beam->intersurface;
  for (k = 1; k < length; k++) {
    beam->rhoL[k] = beam->rho[k-1];
    beam->rhoR[k] = beam->rho[k];
    beam->uL[k] = beam->u[k-1];
    beam->uR[k] = beam->u[k];
    beam->eL[k] = beam->cs[k-1];
    beam->eR[k] = beam->cs[k];
    beam->v_perp_L[0][k] = beam->v_perp[0][k-1];
    beam->v_perp_R[0][k] = beam->v_perp[0][k];
    beam->v_perp_L[1][k] = beam->v_perp[1][k-1];
    beam->v_perp_R[1][k] = beam->v_perp[1][k];
  }
  AdjustBeamBoundaries (beam);
}
