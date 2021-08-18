#include "jupiter.h"

void ViscousStress (dim)	/* This    function    evaluates   the
				   following  three   viscous   stress
				   tensor    components:    Tau_dim_0,
				   Tau_dim_1,   Tau_dim_2, and  stores
				   them  in the StressTensor field  of
				   the current fluid patch */
/* This function evaluates StressTensor[3], which gives the flux of each
   momentum  component   (per    unit  area)  across   the   interface
   perpendicular to dimension dim */
     long dim;
{
  FluidWork *fw;
  long dm[3], i, j, k, m, gncell[3], m_minus;
  long stride[3], iddim, id;
  real *Divergence, *rho, *Radius, *Colatitude, radius=0.0;
  real colatitude = 0.0, centralradius = 0.0;
  real *StressTensor[3], *Edges[3];
  real *InvVolume, *InterSurface[3], *vel[3], *V, *center[3];
  real Dij, Dji;
  real mu;
  real metric[] = {1.0, 1.0, 1.0};
  real Viscosity=0.0;
  fw = CurrentFluidPatch;
  if (fw->desc->level < VISCUTOFFLEVEL)
    Viscosity = VISCOSITY;
  dm[0] = dim;
  dm[1] = (dim == 0);
  dm[2] = 2 - (dim == 2);

		real visc;
	visc = VISCOSITY;//gas
	if (fw->Fluid->next != NULL){
		visc = 0.0; //dust
	}





  Divergence = fw->Divergence;
  InvVolume = fw->desc->InvVolume;
  V = fw->Velocity[dim];
  rho = fw->Density;
  Radius = fw->desc->Center[_RAD_];
  Colatitude = fw->desc->Center[_COLAT_];
  getgridsize (fw->desc, gncell, stride);
  for (j = 0; j < 3; j++) {
    vel[j] = fw->Velocity[j];
    Edges[j] = fw->desc->Edges[j];
    center[j] = fw->desc->Center[j];
    InterSurface[j] = fw->desc->InterSurface[j];
    StressTensor[j] = fw->StressTensor[j];
  }
				/* We first    evaluate   the velocity
				   divergence. The divergence will  be
				   used  once since it features on the
				   stress tensor diagonal only, and we
				   only  evaluate one row or column of
				   that   tensor.   This    divergence
				   therefore needs  to be evaluated at
				   the   center of  the     interfaces
				   perpendicular to dimension `dim` */
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]+(dim == 2); k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]+(dim == 1); j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]+(dim == 0); i++) {
  	m  = i*stride[0]+j*stride[1]+k*stride[2];
	iddim = i;
	if (dim == 1) iddim = j;
	if (dim == 2) iddim = k;
	m_minus = m-stride[dim];
	mu = 2.0*visc*(rho[m]*rho[m_minus])/(rho[m]+rho[m_minus]);
	if (__CYLINDRICAL) {
	  if (dim == _RAD_) radius = Edges[_RAD_][iddim];
	  else radius = Radius[m];
	  metric[_AZIM_] = radius;
	}
	if (__SPHERICAL) {
	  if (dim == _RAD_) centralradius = Edges[_RAD_][iddim];
	  else centralradius = Radius[m];
	  if (dim == _COLAT_) colatitude = Edges[_COLAT_][iddim];
	  else colatitude = Colatitude[m];
	  radius = centralradius * sin(colatitude);
	  metric[_COLAT_] = centralradius;
	  metric[_AZIM_] = radius;
	}
	/* We use the relationship Tau_ij=mu*(Dij+Dji)-2/3*mu*delta_ij*divV +Source Term*/
	/* where Dij = g_i^{-1}g_j\partial_i v_j, where g=1, 1, 1 in cartesian*/
	/* g=1, r, 1 in cylindrical, g=1, r\sin\theta, r in spherical. */
	/* In our notation i is dim, j is id */
	for (id = 0; id < NDIM; id++) {
	  Dij = (vel[id][m]-vel[id][m_minus])/(center[dim][m]-center[dim][m_minus]);
	  Dji =	Dij;
	  if (id != dim) {
	    Dji = .5*((V[m+stride[id]]+V[m_minus+stride[id]])-		\
		      (V[m-stride[id]]+V[m_minus-stride[id]]))/		\
	      (center[id][m+stride[id]]-center[id][m-stride[id]]); /* no average required with m_minus */
	    Dij *= metric[id] / metric[dim];
	    Dji *= metric[dim]/ metric[id];
	  }
	  StressTensor[id][m] = mu*(Dij + Dji -				\
				    (id==dim ? (Divergence[m_minus]+Divergence[m])*ONETHIRD : 0.0));
	  /* Geometry dependent terms are added here */
	  if (__CYLINDRICAL && (dim == _AZIM_) && (id == _AZIM_) && (_RAD_ < NDIM))
	    StressTensor[id][m] += mu*(vel[_RAD_][m]+vel[_RAD_][m_minus])/radius;
	  if (__SPHERICAL && (dim == _AZIM_) && (id == _AZIM_) && (_RAD_ < NDIM))
	    StressTensor[id][m] += mu*(vel[_RAD_][m]+vel[_RAD_][m_minus])/centralradius;
	  if (__SPHERICAL && (dim == _AZIM_) && (id == _AZIM_) && (_COLAT_ < NDIM))
	    StressTensor[id][m] += mu*(vel[_COLAT_][m]+vel[_COLAT_][m_minus])*cos(colatitude)/ \
	      sin(colatitude);
	  if (__SPHERICAL && (dim == _COLAT_) && (id == _COLAT_) && (_RAD_ < NDIM))
	    StressTensor[id][m] += mu*(vel[_RAD_][m]+vel[_RAD_][m_minus])/centralradius;
	}
      }
    }
  }
}
