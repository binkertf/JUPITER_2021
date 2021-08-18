#include "jupiter.h"

void Divergence ()
{
  FluidWork *fw;
  long i, j, k, m, gncell[3], m_minus, m_plus, dim, stride[3];
  real *Divergence, *Radius, *Colatitude;
  real *InvVolume, *InterSurface[3], *vel[3];
  real LowerInterface, UpperInterface, ivol;
  real vplus, vminus, incr_div, radius;
  fw = CurrentFluidPatch;
  Divergence = fw->Divergence;
  InvVolume = fw->desc->InvVolume;
  Radius = fw->desc->Center[_RAD_];
  Colatitude = fw->desc->Center[_COLAT_];
  getgridsize (fw->desc, gncell, stride);
  for (j = 0; j < 3; j++) {
    vel[j] = fw->Velocity[j];
    InterSurface[j] = fw->desc->InterSurface[j];
  }
				/* We     evaluate     the    velocity
				   divergence, centered  at the center
				   of  the zones. For  that we  use an
				   average of  the velocities, and the
				   exact  value of the  zone interface
				   surface area */  
  for (dim = 0; dim < NDIM; dim++) {
    for (k = Nghost[2]-1; k < gncell[2]-Nghost[2]+1; k++) {
      for (j = Nghost[1]-1; j < gncell[1]-Nghost[1]+1; j++) {
	for (i = Nghost[0]-1; i < gncell[0]-Nghost[0]+1; i++) {
	  m  = i*stride[0]+j*stride[1]+k*stride[2];
	  m_minus = m-stride[dim];
	  m_plus = m+stride[dim];
	  ivol = InvVolume[m];
	  LowerInterface   = InterSurface[dim][m];
	  UpperInterface   = InterSurface[dim][m_plus];
	  vplus = vel[dim][m]+vel[dim][m_plus];
	  vminus= vel[dim][m]+vel[dim][m_minus];
	  incr_div = .5*(UpperInterface*vplus-LowerInterface*vminus)*ivol;
	  if (__CYLINDRICAL && (dim == _AZIM_)) {
	    radius = Radius[m];
	    incr_div *= radius;
	  }
	  if (__SPHERICAL && (dim == _COLAT_)) {
	    radius = Radius[m];
	    incr_div *= radius;
	  }
	  if (__SPHERICAL && (dim == _AZIM_)) {
	    radius = Radius[m] * sin(Colatitude[m]);
	    incr_div *= radius;
	  }
	  if (!dim)
	    Divergence[m] = incr_div;
	  else
	    Divergence[m] += incr_div;
	}
      }
    }
  }
}
