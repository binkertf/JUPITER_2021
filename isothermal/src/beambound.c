#include "jupiter.h"

void AdjustBeamBoundaries (beam)
     Beam *beam;
{
  real rhoi, ei, ui, vi=0.0, wi=0.0, x, yg=0.0, zg=0.0;
  real *rhog, *eg, *ug, *vg, *wg, xg;
  long n;
  n = beam->length;
  if (beam->true_bc[INF] > 0) {
    rhoi           = beam->rhoR[NGH];
    ei             = beam->cs[NGH];
    if (Isothermal) ei *= ei;    /* The 'boundary' function below works on cs2, not cs */
    ui             = beam->uR[NGH];
    if (NDIM>1) vi = beam->v_perp_R[0][NGH];
    if (NDIM>2) wi = beam->v_perp_R[1][NGH];
    rhog = &beam->rhoL[NGH];
    eg   = &beam->cs[NGH-1];
    if (Isothermal) *eg *= *eg;  /* Necessary because *eg might be untouched by 'boundary' */
                 /* in which case the 'sqrt' below would result in an */
                 /* erroneous value */
    ug   = &beam->uL[NGH];
    vg   = &beam->v_perp_L[0][NGH];
    wg   = &beam->v_perp_L[1][NGH];
    xg   = beam->rawcoord[NGH-1];
    x    = beam->rawcoord[NGH];
    boundary (rhoi,ei,ui,vi,wi,rhog,eg,\
	      ug,vg,wg,x,xg,yg,zg,beam->true_bc[INF],-1,TRUE);
    if (Isothermal) *eg = sqrt(*eg);
  }
  if (beam->true_bc[SUP] > 0) {
    rhoi           = beam->rhoL[n-NGH];
    ei             = beam->cs[n-NGH-1];
    if (Isothermal) ei *= ei;
    ui             = beam->uL[n-NGH];
    if (NDIM>1) vi = beam->v_perp_L[0][n-NGH];
    if (NDIM>2) wi = beam->v_perp_L[1][n-NGH];
    rhog = &beam->rhoR[n-NGH];
    eg   = &beam->cs[n-NGH];
    if (Isothermal) *eg *= *eg;
    ug   = &beam->uR[n-NGH];
    vg   = &beam->v_perp_R[0][n-NGH];
    wg   = &beam->v_perp_R[1][n-NGH];
    xg   = beam->rawcoord[n-NGH];
    x    = beam->rawcoord[n-NGH-1];
    yg   = beam->rawcoord1;
    zg   = beam->rawcoord2;
    boundary (rhoi,ei,ui,vi,wi,rhog,eg,\
	      ug,vg,wg,x,xg,yg,zg,beam->true_bc[SUP],-1,TRUE);
    if (Isothermal) *eg = sqrt(*eg);
  }
}
