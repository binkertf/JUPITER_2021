#include "jupiter.h"

void FillSources_pot (flag, loc) 
     long flag, loc;
{
  long i, j, k, l, m, smp, smm, Size;
  (void) flag;
  FluidWork *fw;
  tGrid_CPU *d;
  real *inter[3], *center[3], *invvol, *rho, *pot, *v[3];
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
	    if (EXTERNALPOTENTIAL == YES) {
	      smp = m+stride[l];
	      smm = m-stride[l];
	      sinvdx = 1.0/(center[l][smp]-center[l][smm]);
	      sinvdx *= metric_coef;
        if (KEPLERIAN)
          sv[l][m] = pot[m]*pot[m]*(1./pot[smp]-1./pot[smm])*sinvdx; //???
        else
          sv[l][m] = (pot[smm]-pot[smp])*sinvdx;
        } else {
	      sv[l][m] = 0.0;
	    }
	  }
	}
      }
    }
  }
}
