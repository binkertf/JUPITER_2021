#include "jupiter.h"

void PressureCorrection (dt)
     real dt;
{
  FluidWork *fw;
  real *gp[3], *rho;
  real *center[3], *_radius, *_colatitude;
  real *vel[3];
  real *inter[3], *invvol;
  long i, j, k, gncell[3], stride[3], dim;
  long m, mp, l;
  real metric, pp, pm;
  fw = CurrentFluidPatch;
  rho = fw->Density;
  getgridsize (fw->desc, gncell, stride);
  for (dim = 0; dim < 3; dim++) {
    gp[dim] = fw->InterfacePressure[dim];
  }
  for (dim = 0; dim < 3; dim++) {
    center[dim] = fw->desc->Center[dim];
    vel[dim] = fw->Velocity[dim];
    inter[dim] = fw->desc->InterSurface[dim];
  }
  _radius = center[_RAD_];
  _colatitude = center[_COLAT_];
  invvol = fw->desc->InvVolume;
  for (l = 0; l < NDIM; l++) {
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  /* The 'metric' correction is simply because we deal here
	     with the angular velocity in the case specified below*/
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  metric = 1.0;
	  if (__CYLINDRICAL  && (l == _AZIM_))
	    metric = _radius[m];
	  if (__SPHERICAL && (l == _AZIM_))
	    metric = _radius[m]*sin(_colatitude[m]);
	  if (__SPHERICAL && (l == _COLAT_))
	    metric = _radius[m];
	  mp = m+stride[l];
	  pp = gp[l][mp];
	  pm = gp[l][m];
	  vel[l][m] += .5*(pp+pm)*\
	    (inter[l][mp]-inter[l][m])*invvol[m]*dt/(rho[m]*metric);
	  /* The 'metric' correction is due to the fact that we are
	     possibly dealing with an angular velocity here. There is
	     no correction (metric=1.0) in the linear case*/
	}
      }
    }
  }
}

