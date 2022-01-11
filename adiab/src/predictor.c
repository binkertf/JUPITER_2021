#include "jupiter.h"

void Predictor (dt)
     real dt;
{
  long i, j, k, l, m, n, p, Size;
  FluidWork *fw;
  real *rho, *rhop, *e, var;
  real *v[3], *vp[3], *source_v[3], *source_rho;
  real *slope_rho[3], *slope_v[3][3];
  long gncell[3], stride[3], ngh[3];
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho = fw->Density;
  rhop = fw->Density_Pred;
  e = fw->Energy;
  source_rho = fw->SourceRhoPred;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    vp[i] = fw->Velocity_Pred[i];
    source_v[i] = fw->SourceVelocity[i];
    slope_rho[i] = fw->SlopeDensity[i];
    for (j = 0; j < NDIM; j++) {
      slope_v[i][j] = fw->SlopeVelocity[i][j];
    }
  }
  for (i = 0; i < 3; i++)
    ngh[i] = (Nghost[i] > 0 ? Nghost[i]-1 : 0);
  for (l = 0; l < NDIM; l++) {
    for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
      for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	      for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	        m = i*stride[0]+j*stride[1]+k*stride[2];
	        var = 0.0;
	        for (n=0; n < NDIM; n++) {
	          var -= v[n][m]*slope_rho[n][m];
	          var -= rho[m]*slope_v[n][n][m];
	        }
	        var += source_rho[m];
	        rhop[m] = rho[m] + .5*dt*var;
	        for (p=0; p < NDIM; p++) { /* Loop on velocity component */
	          var = 0.0;
	          for (n=0; n < NDIM; n++) /* Loop on slope direction */
	            var += -v[n][m]*slope_v[n][p][m];
	          var += source_v[p][m];
	          var += -e[m]/rho[m]*slope_rho[p][m];
	          /* Note that the -Grad Cs^2 is missing in this prediction step. */
            /* The Cs^2 term is accounted for in FillSources_Predict() and added to source_v */
	          vp[p][m] = v[p][m] + .5*dt*var;
	        }
	      }
      }
    }
  }
}
