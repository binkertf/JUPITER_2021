#include "jupiter.h"

void FillEnergyTot ()
{
  FluidWork *fw;
  long gncell[3], stride[3], i, j, k, size[3];
  long m, idm, len;
  real *Velocity[3], *Density, *InvVolume, *Velf;
  real *energy, *energy_tot;
  real *_radius, *_colatitude;
  real vlin, radius;
  fw = CurrentFluidPatch;
  InvVolume = fw->desc->InvVolume;
  getgridsize (fw->desc, gncell, stride);
  len = 1;
  for (i = 0; i < NDIM; i++) {
    Velocity[i] = fw->Velocity[i];
    size[i] = gncell[i]-2*Nghost[i];
    len *= gncell[i];
  }
  _radius = fw->desc->Center[_RAD_];
  _colatitude = fw->desc->Center[_COLAT_];
  energy = fw->Energy;
  energy_tot = fw->Energy_tot;
  Density = fw->Density;
  memcpy (energy_tot, energy, (size_t)(len*sizeof(real)));
  for (idm = 0; idm < NDIM; idm++) {
    Velf = Velocity[idm];
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  vlin = Velocity[idm][m]; /* Linear velocity *in corotating frame* */
	  if (__CYLINDRICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    vlin = radius*vlin;
	  }
	  if (__SPHERICAL && (idm == _AZIM_)) {
	    radius = _radius[m];
	    radius *= sin(_colatitude[m]);
	    vlin = radius*vlin;
	  }
	  if (__SPHERICAL && (idm == _COLAT_)) { /* The colatitude velocity is angular */
	    radius = _radius[m];
	    vlin *= radius;
	  }
	  energy_tot[m] += Density[m]*.5*vlin*vlin;	/* Volumic total energy */
	}
      }
    }
  }
}
