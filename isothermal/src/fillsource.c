#include "jupiter.h"

void FillSources (flag, loc)
     long flag, loc;
{
  long i, j, k, l, m, mp[3], mm[3], smp, smm, Size;
  FluidWork *fw;
  tGrid_CPU *d;
  real *inter[3], *center[3], *invvol, *rho, *pot, *v[3];
  real *sv[3], *met[3][2], *acc[3];
  real metric_coef=1.0;
  real sinvdx, vphi, vtheta, rad, cot, sine, *invm[3][2];
  real *srcrho, dsurf, *sinecol;
  long gncell[3], stride[3], ngh[3], ind[3], ip1, ip2;
  fw = CurrentFluidPatch;
  d = fw->desc;
  getgridsize (d, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho = fw->Density;
  pot = fw->Potential;
  invvol = d->InvVolume;
  srcrho = fw->SourceRhoPred;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    sv[i] = fw->SourceVelocity[i];
    acc[i] = fw->Accel[i];
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
  if (flag == PREDICT) {
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
		  sv[l][m] = pot[m]*pot[m]*(1./pot[smp]-1./pot[smm])*sinvdx;
		else
		  sv[l][m] = (pot[smm]-pot[smp])*sinvdx;
	      } else {
		sv[l][m] = 0.0;
	      }
	      if (GridFriction[l] != 0.0)
		sv[l][m] -= v[l][m]*GridFriction[l];
	    }
	  }
	}
      }
    }
  }

  if (flag == PREDICT)
    memcpy (acc[0], sv[0], NDIM*Size*sizeof(real));
  else
    memcpy (sv[0], acc[0], NDIM*Size*sizeof(real));

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
	      if (_COLAT_ < NDIM) //azimuthal direction
		      vtheta = v[_COLAT_][m];

	      if (_RAD_ < NDIM) //radial direction
		      sv[_RAD_][m] += rad*(vphi*vphi+vtheta*vtheta); // source velcoity

	      if (_COLAT_ < NDIM) { //polar direction
		        cot = (inter[_COLAT_][mp[_COLAT_]]-inter[_COLAT_][m])*invvol[m]*rad;
		        sv[_COLAT_][m] += rad*vphi*vphi*cot;
	          }

	    } else { //not spherical
	      if (_RAD_ < NDIM)
		    sv[_RAD_][m] += rad*vphi*vphi;
	    }
	  }
	}
      }
    }
  }
  if (flag == PREDICT) {
    for (i = 0; i < Size; i++)
      srcrho[i] = 0.0;
    for (l = 0; l < NDIM; l++) {
      ip1 = (l == 0);
      ip2 = 2 - (l == 2);
      for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
	ind[2] = k;
	for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	  ind[1] = j;
	  for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	    ind[0] = i;
	    m = i*stride[0]+j*stride[1]+k*stride[2];
	    if (d->CommSource[m] & loc) {
	      if (!__CARTESIAN) {
		smp = m+stride[l];
		metric_coef = met[l][0][ind[ip1]]*met[l][1][ind[ip2]];
		dsurf = inter[l][smp]-inter[l][m];
		srcrho[m] += v[l][m]*dsurf*metric_coef;
	      }
	    }
	  }
	}
      }
    }
    for (i = 0; i < Size; i++)
      srcrho[i] *= rho[i]*invvol[i];
  }
}



void FillDiffSources_Dust (flag, loc)
     long flag, loc;
{
  long i, j, k, l, m, mp[3], mm[3], smp, smm, Size;
  FluidWork *fw;
  tGrid_CPU *d;
  real *inter[3], *center[3], *invvol, *rho, *rho_g, *cs2, *pot, *v[3], *a2;
  real *sv[3], *met[3][2], *acc[3];
  real metric_coef=1.0;
  real delta, sinvdx, vphi, vtheta, rad, cot, sine, *invm[3][2];
  real *srcrho, dsurf, *sinecol;
  long gncell[3], stride[3], ngh[3], ind[3], ip1, ip2;
  fw = CurrentFluidPatch;
  d = fw->desc;
  getgridsize (d, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho = fw->Density;
  a2 = fw->Energy; //diffusion sound speed squared

  rho_g = fw->next->Density->Field;
  cs2 = fw->next->Energy->Field;

  pot = fw->Potential;
  invvol = d->InvVolume;
  srcrho = fw->SourceRhoPred;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    sv[i] = fw->SourceVelocity[i];
    acc[i] = fw->Accel[i];
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
  if (flag == PREDICT) {
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

    sv[l][m] +=  1. / rho_g[m] * (a2[smp] * rho_g[smp] - a2[smm] * rho_g[smm])*sinvdx;



	    }
	  }
	}
      }
    }
  }

}





void FillDiffSources_Gas (flag, loc)
     long flag, loc;
{
  long i, j, k, l, m, mp[3], mm[3], smp, smm, Size;
  FluidWork *fw;
  tGrid_CPU *d;
  real *inter[3], *center[3], *invvol, *rho_d, *rho_g, *cs2, *pot, *v[3], *a2;
  real *sv[3], *met[3][2], *acc[3];
  real metric_coef=1.0;
  real delta, sinvdx, vphi, vtheta, rad, cot, sine, *invm[3][2];
  real *srcrho, dsurf, *sinecol;
  long gncell[3], stride[3], ngh[3], ind[3], ip1, ip2;
  fw = CurrentFluidPatch;
  d = fw->desc;
  getgridsize (d, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho_g = fw->Density;
  cs2 = fw->Energy; //gas sound speed squared

  rho_d = fw->prev->Density->Field;
  a2 = fw->prev->Energy->Field;

  pot = fw->Potential;
  invvol = d->InvVolume;
  srcrho = fw->SourceRhoPred;
  for (i = 0; i < NDIM; i++) {
    v[i] = fw->Velocity[i];
    sv[i] = fw->SourceVelocity[i];
    acc[i] = fw->Accel[i];
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
  if (flag == PREDICT) {
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

    sv[l][m] -= a2[m] * (rho_d[smp] / rho_g[smp] - rho_d[smm] /rho_g[smm])*sinvdx;

	    }
	  }
	}
      }
    }
  }

}
