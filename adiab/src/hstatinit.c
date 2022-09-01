/** \file hstatinit.c
    Tabulate the values of the equilibrium quantities.
    Tabulate the values of the zone centered density, and of the density,
    sound speed and pressure at the interface, for the equilibrium that
    one aims at enforcing. The interface centered sound speed and pressure
    are obtained through an integration on the interface, using a trapezium
    method and a number of elements that depend on the resolution, so that
    a value on a coarse interface coincides with the properly averaged values
    of the corresponding fine interfaces. Note that by doing so we do not
    have exactly the relationship P=rho*cs^2 at the interface, but this is
    not a problem. Note that rho at the interface is obtained by direct
    evaluation of the analytic expression of the equilibrium state
    (hence the inclusion of libinit.cx)

*/

#include "jupiter.h"
#include "init.h"

void SetCenteredDensity (fluid)
     FluidPatch *fluid;
{
  SqueezedField *rhoc, *enec;
  long gncell[3], i, j, k, ms, m;
  long ii, jj, kk;
  long strides[3], stride[3];
  rhoc = fluid->Rho_eq_c;
  enec = fluid->Ene_eq_c;
  for (i = 0; i < 3; i++) {
    gncell[i] = rhoc->gncell[i];
    stride[i] = fluid->desc->stride[i];
    strides[i]= rhoc->stride[i];
  }
  /* Note that 'gncell' below refers to that of the squeezed array */
  for (k = 0; k < gncell[2]; k++) {
    for (j = 0; j < gncell[1]; j++) {
      for (i =0; i < gncell[0]; i++) {
	/* The shift below is meant to use values in the active part
	   of the mesh along the squeezed dimensions */
	ii = i; jj = j; kk = k;
	if (SqueezeDim[0]) ii+= Nghost[0];
	if (SqueezeDim[1]) jj+= Nghost[1];
	if (SqueezeDim[2]) kk+= Nghost[2];
	ms = ii*strides[0]+jj*strides[1]+kk*strides[2];
	m  = ii*stride[0]+jj*stride[1]+kk*stride[2];
	rhoc->field[ms] = fluid->Density->Field[m];
	if (!Isothermal)
	  enec->field[ms] = fluid->Energy->Field[m];
      }
    }
  }  
}

void SetInterfaceQuantities (fluid)
/* Here we set the density and pressure at the interface, that
   correspond to a hydrostatic equilibrium */
     FluidPatch *fluid;
{
  SqueezedField *rhoi[3], *pi[3], *cs2i[3], *enei[3];
  long gncell[3], i, j, k, ms, m, l, idx[3], nsteps;
  long strides[3], stride[3];
  real corner_min[3], corner_max[3];
  real radius, azimuth, colatitude, x, y, z;
  real vx, vy, vz, vrad, v_azimuth, v_colatitude;
  real density, energy, axial_radius, gamma, pressure;
  real interm1, interm2, interm3; /* Work variables for user convenience */
  long *ml, n1, n2, mc, mlev1, mlev2;
  SetFluidProperties (fluid);
  for (l = 0; l < 3; l++) {
    rhoi[l] = fluid->Rho_eq_i[l];
    pi[l] = fluid->Pres_eq_i[l];
    cs2i[l] = fluid->Cs2_i[l];
    if (!Isothermal)
      enei[l] = fluid->Ene_eq_i[l];
    stride[l] = fluid->desc->stride[l];
    if (CorrHydroStat[l]) {
      for (i = 0; i < 3; i++) {
	gncell[i] = rhoi[l]->gncell[i];
	strides[i]= rhoi[l]->stride[i];
      }
    }
  }
  if (KEPLERIAN) {
    ml = fluid->desc->MaxLevMeridianProj;
    n1 = gncell[_RAD_];
    n2 = gncell[_COLAT_];
  }
  for (l = 0; l < 3; l++) {
    if (CorrHydroStat[l]) {
      for (k = 0; k < gncell[2]; k++) {
		idx[2] = k;
		for (j = 0; j < gncell[1]; j++) {
	  	  	idx[1] = j;
	  	  	for (i = 0; i < gncell[0]; i++) {
	    		idx[0] = i;
	    		ms = i*strides[0]+j*strides[1]+k*strides[2];
	    		m = i*stride[0]+j*stride[1]+k*stride[2];
	    		/* We need to know the physical coordinates of the
	       		interface in order to evaluate the density at the
	       		latter */
	    		radius = x = fluid->desc->Center[_RAD_][m];
	    		if (_RAD_ == l)
	      			radius = x = fluid->desc->Edges[_RAD_][idx[_RAD_]];
	    		azimuth = y = fluid->desc->Center[_AZIM_][m];
	    		if (_AZIM_ == l) 
	      			azimuth = y = fluid->desc->Edges[_AZIM_][idx[_AZIM_]];
	    		colatitude = z = fluid->desc->Center[_COLAT_][m];
	    		if (_COLAT_ == l)
	      			colatitude = z = fluid->desc->Edges[_COLAT_][idx[_COLAT_]];
	    		density = energy = -1e20;

				#include "libinit.cx"

				if (!Isothermal && Stellar && !CONST_GAMMA){
      				pressure = secant_method(energy, density, 0.00001); // pressure
      				gamma = Gamma1(pressure * MU / density * TEMP0, density);
    			}
				else
					gamma = GAMMA;

	    		rhoi[l]->field[ms] = density;
	    		if (!Isothermal)
	      			enei[l]->field[ms] = energy;
	    		corner_min[0] = fluid->desc->Edges[_RAD_][idx[_RAD_]];
	    		corner_min[1] = fluid->desc->Edges[_AZIM_][idx[_AZIM_]];
	    		corner_min[2] = fluid->desc->Edges[_COLAT_][idx[_COLAT_]];
	    		corner_max[0] = fluid->desc->Edges[_RAD_][idx[_RAD_]+1];
	    		corner_max[1] = fluid->desc->Edges[_AZIM_][idx[_AZIM_]+1];
	    		corner_max[2] = fluid->desc->Edges[_COLAT_][idx[_COLAT_]+1];
	    		nsteps = 2L<<(LevMax-fluid->desc->level);
	    		if (KEPLERIAN) {
	      			mc = idx[_RAD_]+idx[_COLAT_]*n1;
	      			if (_COLAT_ < _RAD_)
						mc = idx[_COLAT_]+idx[_RAD_]*n2;
	      			if (l == _COLAT_) {
						mlev1 = ml[mc];
						if (idx[_COLAT_] > 0) {
		  					mlev2=ml[mc-(_RAD_ < _COLAT_ ?  n1 : 1)];
		  					if (mlev2 > mlev1) mlev1 = mlev2;
						}	 
	      			}
	    			if (l == _RAD_) {
						mlev1 = ml[mc];
						if (idx[_RAD_] > 0) {
		  					mlev2=ml[mc-(_RAD_ < _COLAT_ ?  1 : n2)];
		  					if (mlev2 > mlev1) mlev1 = mlev2;
						}	 
	      			}
	      			nsteps = 2L<<(1+mlev1-fluid->desc->level);
	    		}
	    		if (NDIM == 3) {
	      			if (CoordNb[l] != 2) {
						pi[l]->field[ms] = IC_2D_Mean (corner_min, corner_max, IC_Pressure, nsteps, CoordNb[l]);
	      			}
	      			else {		/* In colatitude, we use a 'radius' weighting in the average */
						pi[l]->field[ms] = IC_2D_Mean (corner_min, corner_max, IC_PressureByRad, nsteps, 2);
						pi[l]->field[ms] /= fluid->desc->Center[_RAD_][m];
	      			}
	    		}
	    		if (NDIM == 2) {
	      			if (l == _RAD_)   corner_min[0] = corner_max[0] = radius;
	      			if (l == _AZIM_)  corner_min[1] = corner_max[1] = azimuth;
	      			if (l == _COLAT_) corner_min[2] = corner_max[2] = colatitude;
	      			pi[l]->field[ms] = IC_Average (corner_min, corner_max, IC_Pressure, nsteps);
	    		}
	    		if (NDIM == 1) {
	      			if (Isothermal)
						pi[l]->field[ms] = density * energy;
	      			else
						pi[l]->field[ms] = energy*(gamma-1.);
	    		}
	    		cs2i[l]->field[ms] = pi[l]->field[ms]/density;
	    		if (!Isothermal)
	      			cs2i[l]->field[ms] *= gamma;
	  		}
		}
      }   
    }
  }
}
