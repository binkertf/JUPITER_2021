#include "jupiter.h"


void DiskRadialDrift(dt)
     real dt;
{
  FluidWork *fw;
  long i, j, k, stride[3], gncell[3], m;
  tGrid_CPU *desc;
  real radius;
  real *vazimuth;
  
  fw = CurrentFluidPatch;
  desc = fw->desc;
  getgridsize (desc, gncell, stride);
  vazimuth = fw->Velocity[_AZIM_];
  
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	radius = desc->Center[_RAD_][m];
	vazimuth[m] += dt * DRIFTVELOCITY * 0.5 * pow(radius,-(3.5+SIGMASLOPE));
      }
    }
  }
}
