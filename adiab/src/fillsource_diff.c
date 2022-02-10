#include "jupiter.h"

void FillSources_diff_dust (flag, loc, dt) 
     long flag, loc;
     real dt;
{
  long i, j, k, l, m, smp, smm, Size;
  (void) flag;
  FluidWork *fw;
  tGrid_CPU *d;
  real *inter[3], *center[3], *invvol, *rho, *pot, *v[3], *a2, *rho_g, *cs2;
  real *sv[3], *met[3][2];
  real metric_coef=1.0;
  real sinvdx, *invm[3][2];
  real *sinecol;
  long gncell[3], stride[3], ngh[3];
  fw = CurrentFluidPatch;
  d = fw->desc;
  getgridsize (d, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho = fw->Density;
  pot = fw->Potential;
  a2 = fw->Energy; //diffusion sound speed squared

  rho_g = fw->next->Density->Field;
  cs2 = fw->next->Energy->Field; //gas sound speed squared (isothermal) / internal energy (radiative)

  invvol = d->InvVolume;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    sv[i] = fw->SourceVelocity[i];
  }
  for (i = 0; i < 3; i++) {
    inter[i] = d->InterSurface[i];
    center[i] = d->Center[i];
    invm[i][0] = d->InvMetric[i][0];
    invm[i][1] = d->InvMetric[i][1];
    met[i][0] = d->Metric[i][0];
    met[i][1] = d->Metric[i][1];
  }
  sinecol = met[_AZIM_][_COLAT_ > _RAD_];
  for (i = 0; i < 3; i++)
    ngh[i] = (Nghost[i] > 0 ? Nghost[i]-1 : 0);
  for (l = 0; l < NDIM; l++) {
    for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
      for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	      if (l == 0) metric_coef = invm[0][0][j] * invm[0][1][k];
	      for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	        if (l == 2) metric_coef = invm[2][0][i] * invm[2][1][j];
	        if (l == 1) metric_coef = invm[1][0][i] * invm[1][1][k];
	        m = i*stride[0]+j*stride[1]+k*stride[2];
	        if (d->CommSource[m] & loc) {
	          smp = m+stride[l];
	          smm = m-stride[l];
	          sinvdx = 1.0/(center[l][smp]-center[l][smm]);
	          sinvdx *= metric_coef;
	          sv[l][m] = 1. / rho_g[m] * (a2[smp] * rho_g[smp] - a2[smm] * rho_g[smm]) * sinvdx;

            //diffusion limiter
            if ((sv[l][m]*dt)>a2[m]){
              sv[l][m] = a2[m] / dt; 
            } 
            if ((sv[l][m]*dt)<-a2[m]){
              sv[l][m] = -a2[m] / dt; 
            }
	        }
	      }
      }
    }
  }
}


void FillSources_diff_gas (flag, loc, dt) 
     long flag, loc;
     real dt;
{
  long i, j, k, l, m, smp, smm, Size;
  (void) flag;
  FluidWork *fw;
  tGrid_CPU *d;
  real *inter[3], *center[3], *invvol, *rho_d, *pot, *v[3], *a2, *rho_g, *cs2;
  real *sv[3], *met[3][2];
  real metric_coef=1.0;
  real sinvdx, *invm[3][2];
  real *sinecol;
  long gncell[3], stride[3], ngh[3];
  fw = CurrentFluidPatch;
  d = fw->desc;
  getgridsize (d, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho_d = fw->Fluid->prev->Density->Field; //dust density
  a2 = fw->Fluid->prev->Energy->Field; //diffusion sound speed squared
  rho_g = fw->Density; //gas density
  cs2 = fw->Energy; //gas sound speed squared (isothermal) / internal energy (radiative)

  invvol = d->InvVolume;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    sv[i] = fw->SourceVelocity[i];
  }
  for (i = 0; i < 3; i++) {
    inter[i] = d->InterSurface[i];
    center[i] = d->Center[i];
    invm[i][0] = d->InvMetric[i][0];
    invm[i][1] = d->InvMetric[i][1];
    met[i][0] = d->Metric[i][0];
    met[i][1] = d->Metric[i][1];
  }
  sinecol = met[_AZIM_][_COLAT_ > _RAD_];
  for (i = 0; i < 3; i++)
    ngh[i] = (Nghost[i] > 0 ? Nghost[i]-1 : 0);
  for (l = 0; l < NDIM; l++) {
    for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
      for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	      if (l == 0) metric_coef = invm[0][0][j] * invm[0][1][k];
	      for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	        if (l == 2) metric_coef = invm[2][0][i] * invm[2][1][j];
	        if (l == 1) metric_coef = invm[1][0][i] * invm[1][1][k];
	        m = i*stride[0]+j*stride[1]+k*stride[2];
	        if (d->CommSource[m] & loc) {
	          smp = m+stride[l];
	          smm = m-stride[l];
	          sinvdx = 1.0/(center[l][smp]-center[l][smm]);
	          sinvdx *= metric_coef;
	          sv[l][m] = a2[m] * (rho_d[smp] / rho_g[smp] - rho_d[smm] / rho_g[smm]) * sinvdx;
            
            //diffusion limiter
            if ((sv[l][m]*dt)>a2[m]){
              sv[l][m] = a2[m] / dt; 
            } 
            if ((sv[l][m]*dt)<-a2[m]){
              sv[l][m] = -a2[m] / dt; 
            }
	        }
	      }
      }
    }
  }
}
