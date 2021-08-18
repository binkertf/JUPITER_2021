#include "jupiter.h"

void EnergyCorrection (dt)
     real dt;
{
  /* Here we subtract the kinetic energy and power of real forces from
     the total energy. */
  long i,j,k,l,m, gncell[3], stride[3], dim;
  real vlin, temperature, *energy;
  real *radius, *colatitude;
  real *energy_tot, *density, *vel[3],*azim;
  FluidWork *fw;
  (void) dt;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  for (dim = 0; dim < 3; dim++) {
    vel[dim] = fw->Velocity[dim];
  }
  energy = fw->Energy;
  energy_tot = fw->Energy_tot;
  density = fw->Density;
  radius = fw->desc->Center[_RAD_];
  colatitude = fw->desc->Center[_COLAT_];
  azim = fw->desc->Center[_AZIM_];
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	for (l = 0; l < NDIM; l++) {
	  vlin = vel[l][m];
	  if ((l == _AZIM_) && (__CYLINDRICAL))
	    vlin = vlin*radius[m];
	  if ((l == _AZIM_) && (__SPHERICAL))
	    vlin = vlin*radius[m]*sin(colatitude[m]);
	  if ((l == _COLAT_) && (__SPHERICAL))
	    vlin = vlin*radius[m];
/*     if  (fw->desc->level == LevMax &&azim[m] < 2*0.0059107 &&azim[m] > -2*0.0059107 && radius[m] < 1.0+2*0.0059107 && radius[m] > 1.0-2*0.0059107 && colatitude[m] > 1.5704-2*0.0059107 && colatitude[m] < 1.5708) printf("%lg %lg %lg %lg %lg\n", energy_tot[m],.5*density[m]*vlin*vlin,radius[m],colatitude[m],azim[m]);*/
	  energy_tot[m] -= .5*density[m]*vlin*vlin;
	}
#ifdef TEMPFLOOR
	temperature = energy_tot[m]/(CV*density[m]);
	if (temperature > BOUNDTEMP/TEMP0/TEMPRATIO)
#endif
	  energy[m] = energy_tot[m];
#ifdef TEMPFLOOR
	else
	  energy[m] = BOUNDTEMP/TEMP0/TEMPRATIO * CV * density[m];
#endif
      }
    }
  }
}


void EnergyCorrection2 (dt)
     real dt;
{
  /* Here we subtract the kinetic energy and power of real forces from
     the total energy. */
  long i,j,k,l,m, gncell[3], stride[3], dim;
  real vlin, *energy, ur, temperature, metric;
  real *radius, *colatitude, *acc[3];
  real *energy_tot, *density, *vel[3],*azim;
  FluidWork *fw;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  for (dim = 0; dim < 3; dim++) {
    vel[dim] = fw->Velocity[dim];
    acc[dim] = fw->SourceVelocity[dim];
  }
  energy = fw->Energy;
  energy_tot = fw->Energy_tot;
  density = fw->Density;
  radius = fw->desc->Center[_RAD_];
  colatitude = fw->desc->Center[_COLAT_];
  azim = fw->desc->Center[_AZIM_];
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	ur = 0.0;
	for (l = 0; l < NDIM; l++) {
	  vlin = vel[l][m];
	  if ((l == _AZIM_) && (__CYLINDRICAL))
	    vlin = vlin*radius[m];
	  if ((l == _AZIM_) && (__SPHERICAL))
	    vlin = vlin*radius[m]*sin(colatitude[m]);
	  if ((l == _COLAT_) && (__SPHERICAL)) {
	    vlin = vlin*radius[m];
	    ur += vlin*cos(colatitude[m]);
	  }
	  if ((l == _RAD_) && (__SPHERICAL))
	    ur += vlin*sin(colatitude[m]); /* upon loop completion ur
					      is the cylindrical
					      radial velocity */
 /*     if  (fw->desc->level == LevMax &&azim[m] < 2*0.0059107 &&azim[m] > -2*0.0059107 && radius[m] < 1.0+2*0.0059107 && radius[m] > 1.0-2*0.0059107 && colatitude[m] > 1.5704-2*0.0059107 && colatitude[m] < 1.5708) printf("ene1=%lg kin_e=%lg rad=%lg co=%lg az=%lg\n", energy_tot[m],.5*density[m]*vlin*vlin,radius[m],colatitude[m],azim[m]);*/
  /*   if  (fw->desc->level == LevMax &&azim[m] < 2*0.0059107 &&azim[m] > -2*0.0059107 && radius[m] < 1.0+2*0.0059107 && radius[m] > 1.0-2*0.0059107 && colatitude[m] > 1.5704-2*0.0059107 && colatitude[m] < 1.5708) printf("%lg %lg %lg %lg %lg\n", energy_tot[m],.5*density[m]*vlin*vlin,radius[m],colatitude[m],azim[m]);*/
	  energy_tot[m] -= .5*density[m]*vlin*vlin;
	  energy_tot[m] += density[m]*acc[l][m]*vlin*dt; /* Work of gravitational forces */
	}
	/* Work of effective potential (centrifugal term) */
	energy_tot[m] += density[m]*OMEGAFRAME*OMEGAFRAME*sin(colatitude[m])*radius[m]*ur*dt;
#ifdef TEMPFLOOR
	temperature = energy_tot[m]/(CV*density[m]);
	if (temperature > BOUNDTEMP/TEMP0/TEMPRATIO)
#endif
	  energy[m] = energy_tot[m];
#ifdef TEMPFLOOR
	else
	  energy[m] = BOUNDTEMP/TEMP0/TEMPRATIO * CV * density[m];
#endif
      }
    }
  }
}

