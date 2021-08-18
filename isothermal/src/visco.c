#include "jupiter.h"

void ApplyViscousStress (dt)
     real dt;
{
  long dim, ii, i[3], m, m_other_side, gncell[3], stride[3], mdim, ex[3];
  real *momentum[3];
  real radius_inf = 0.0, radius_sup = 0.0, *st[3];
  real centralradius_inf = 0.0, centralradius_sup = 0.0;
  real colatitude_inf = 0.0, colatitude_sup = 0.0;
  real InterfaceInf, lower_face_flux;
  real *Radius, *Colatitude;
  FluidWork *fp;
  fp = CurrentFluidPatch;
  getgridsize (fp->desc, gncell, stride);
  for (ii = 0; ii < 3; ii++) {
    momentum[ii] = fp->Momentum[ii];
    st[ii] = fp->StressTensor[ii];
  }
  Radius = fp->desc->Center[_RAD_];
  Colatitude = fp->desc->Center[_COLAT_];
  Divergence ();
  for (dim = 0; dim < NDIM; dim++) {
    ViscousStress (dim);
    ex[0] = ex[1] = ex[2] = 0;
    ex[dim]=1;
    for (i[2] = Nghost[2]; i[2] < gncell[2]-Nghost[2]+ex[2]; i[2]++) {
      for (i[1] = Nghost[1]; i[1] < gncell[1]-Nghost[1]+ex[1]; i[1]++) {
	for (i[0] = Nghost[0]; i[0] < gncell[0]-Nghost[0]+ex[0]; i[0]++) {
	  m = i[0]*stride[0]+i[1]*stride[1]+i[2]*stride[2];
	  m_other_side = m+stride[dim];
	  InterfaceInf = fp->desc->InterSurface[dim][m];
	  if (__CYLINDRICAL) {
	    if (dim == _RAD_) {
	      radius_inf = fp->desc->Edges[dim][i[dim]];
	      radius_sup = fp->desc->Edges[dim][i[dim]+1];
	    } else {
	      radius_inf = radius_sup = Radius[m];
	    }
	  }
	  if (__SPHERICAL) {
	    if (dim == _RAD_) {
	      centralradius_inf = fp->desc->Edges[dim][i[dim]];
	      centralradius_sup = fp->desc->Edges[dim][i[dim]+1];
	    } else {
	      centralradius_inf = centralradius_sup = Radius[m];
	    }
	    if (dim == _COLAT_) {
	      colatitude_inf = fp->desc->Edges[dim][i[dim]];
	      colatitude_sup = fp->desc->Edges[dim][i[dim]+1];
	    } else {
	      colatitude_inf = colatitude_sup = Colatitude[m];
	    }
	    radius_inf = centralradius_inf * sin(colatitude_inf);
	    radius_sup = centralradius_sup * sin(colatitude_sup);
	  }
	  for (mdim = 0; mdim < NDIM; mdim++) {
	    lower_face_flux = st[mdim][m]*InterfaceInf;
	    if (__CYLINDRICAL && (mdim == _AZIM_)) { /* We work on the angular momentum */
	      lower_face_flux *= radius_inf;
	    }
	    if (__SPHERICAL && (mdim == _AZIM_)) { /* We work on the angular momentum */
	      lower_face_flux *= radius_inf;
	    }
	    if (__SPHERICAL && (mdim == _COLAT_)) { /* We work on r^2\omega_\theta */
	      lower_face_flux *= centralradius_inf;
	    }
//if(fp->desc->level==7&&m==59354&&CPU_Rank==3&&mdim==1) printf("VISCO #1 flux mom=%lg dimm=%ld \n",fp->Flux_mom[1][dim][m],dim);
	    fp->Flux_mom[mdim][dim][m] -= lower_face_flux*dt;
//if(fp->desc->level==7&&m==59354&&CPU_Rank==3&&mdim==1) printf("VISCO #2 flux mom=%lg dimm=%ld st=%lg interface=%lg dt=%lg\n",fp->Flux_mom[1][dim][m],dim,st[1][m],InterfaceInf,dt);
	  }
	  /* Non-conservative curvature terms in curvilinear coordinates are
	     applied below */
	  if (__CYLINDRICAL) {
	    if (dim == _AZIM_) {
	      if (_RAD_ < NDIM)
		momentum[_RAD_][m] -= dt*.5*(st[_AZIM_][m]+\
						st[_AZIM_][m_other_side])/radius_sup;
	      /* radius_sup == radius_inf in the circumstances where the above line
		 applies, so no averaging is required */
	    }
	  }
	  if (__SPHERICAL) {
	    if (dim == _AZIM_) {
	      if (_RAD_ < NDIM)
		momentum[_RAD_][m] -= dt*.5*(st[_AZIM_][m]+\
						st[_AZIM_][m_other_side])/centralradius_sup;
	      /* centralradius_sup == centralradius_inf in the circumstances where the above line
		 applies, so no averaging is required */
	      if (_COLAT_ < NDIM)
		momentum[_COLAT_][m] -= dt*.5*(st[_AZIM_][m]+\
						    st[_AZIM_][m_other_side])*\
		  cos(colatitude_inf)/radius_sup;
	      /* radius_sup == radius_inf in the circumstances where the above line
		 applies, so no averaging is required. Same for colatitude. */
	    }
	    if (dim == _COLAT_) {
	      if (_RAD_ < NDIM)
		momentum[_RAD_][m] -= dt*.5*(st[_COLAT_][m]+\
						st[_COLAT_][m_other_side])/centralradius_sup;
	      /* centralradius_sup == centralradius_inf in the circumstances where the above line
		 applies, so no averaging is required */
	    }
	  }
	}
      }
    }
  }
}
