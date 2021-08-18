#include "jupiter.h"

void FillSources_geom_col (flag, loc) 
     long flag, loc;
{
  long i, j, k, l, m, mp[3], mm[3], Size;
  FluidWork *fw;
  tGrid_CPU *d;
  (void) flag;
  real *inter[3], *center[3], *invvol, *rho, *pot, *v[3];
  real *sv[3], *met[3][2];
  real vphi, vtheta, rad, cot, sine, *invm[3][2];
  real *sinecol;
  long gncell[3], stride[3], ngh[3], ind[3];
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
	for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  sv[l][m] = 0.0;
	}
      }
    }
  }
  for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
    ind[2] = k;
    for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
      ind[1] = j;
      for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	ind[0] = i;
	m = i*stride[0]+j*stride[1]+k*stride[2];
	if (d->CommSource[m] & loc) {
	  for (l = 0; l < 3; l++) { /* 3, not NDIM */
	    mp[l] = m+stride[l];
	    mm[l] = m-stride[l];
	  }
	  if (!__CARTESIAN) {
	    vphi = vtheta = 0.0;
	    if (_AZIM_ < NDIM)
	      vphi = v[_AZIM_][m];
	    vphi += OMEGAFRAME;
	    rad = center[_RAD_][m];
	    if (__SPHERICAL) {
	      sine = sinecol[ind[_COLAT_]];
	      vphi *= sine;
	      if (_COLAT_ < NDIM)
		vtheta = v[_COLAT_][m];
	      if (_COLAT_ < NDIM) {
		cot = (inter[_COLAT_][mp[_COLAT_]]-inter[_COLAT_][m])*invvol[m]*rad;	
		sv[_COLAT_][m] += rad*vphi*vphi*cot;
	      }
	    }
	  }
	}
      }
    }
  }
}
