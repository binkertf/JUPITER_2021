#include "jupiter.h"
#include "init.h"

#define PERCENTMARGIN 30
#define HIGHPERCENTMARGIN 2
/* The above means that we damp perturbations up to +/- 20% of the
   mesh radial boundaries */

void ApplyStockholmBoundaryConditions (dt)
real dt;
{
  FluidWork *fw;
  real *center[3];
  real *vel[3], *dens, ramp=0.0, lambda, R_in, R_out;
  real Tau=0.0, rho0, vrad0, vtheta0, Rmin, Rmax;
  real margin, highmargin, Colmin, Colmax, damp_zone_min, damp_zone_max;
  long i, j, k, m, gncell[3], stride[3], dim;
  boolean correct;
  real radius, azimuth, colatitude, x, y, z;
  real vx, vy, vz, vrad, v_azimuth, v_colatitude;
  real Bx, By, Bz, Brad, B_azimuth, B_colatitude;
  real density, energy, axial_radius, a2m, a2p;
  real interm1, interm2, interm3; /* Work variables for user convenience */
  margin = 0.01*PERCENTMARGIN;
  highmargin = 0.01*HIGHPERCENTMARGIN;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  SetFluidProperties (fw->Fluid);
  for (dim = 0; dim < 3; dim++) {
    center[dim] = fw->desc->Center[dim];
    vel[dim] = fw->Velocity[dim];
  }
  Rmin = corner_min0[_RAD_];
  Rmax = corner_max0[_RAD_];
  Colmin = corner_min0[_COLAT_];
  Colmax = corner_max0[_COLAT_];
  damp_zone_min = margin*fabs(Colmin-M_PI/2.0);
  if (ASPECTRATIO < damp_zone_min) damp_zone_min = ASPECTRATIO;
  damp_zone_max=margin*fabs(Colmax-M_PI/2.0);
  if (ASPECTRATIO < damp_zone_max) damp_zone_max = ASPECTRATIO;
  R_in = Rmin*(1.+margin);
  R_out= Rmax/(1.+margin);
  /*damp_zone_min_high = highmargin*fabs(Colmin-M_PI/2.0);
  damp_zone_max_high = highmargin*fabs(Colmax-M_PI/2.0);
  R_in_high = Rmin*(1.+marginhigh);
  R_out_high= Rmax/(1.+marginhigh);*/
  dens = fw->Density;
  /*  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
    for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {*/
  for (k = 0; k < gncell[2]; k++) {
    for (j = 0; j < gncell[1]; j++) {
      for (i = 0; i < gncell[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	radius = center[_RAD_][m];
	colatitude = center[_COLAT_][m];
	correct = NO;
	x = y = 0.0;
	if (radius < R_in) {
	  correct = YES;
	  x = (R_in-radius)/(R_in-Rmin);
	}
	if (radius > R_out) {
	  correct = YES;
	  x = (radius-R_out)/(Rmax-R_out);
	}
	if (__SPHERICAL && (_COLAT_ < NDIM)) {
	  if (fabs(Colmin - colatitude) < damp_zone_min) {
	    correct = YES;
	    y = 1. - (colatitude - Colmin) / damp_zone_min;
	  }
	  if (fabs(Colmax - colatitude) < damp_zone_max) {
	    correct = YES;
	    y = 1. - (Colmax - colatitude) / damp_zone_max;
	  }
	}
	if (correct == TRUE) {
	  ramp = 1.-(1.-x)*(1.-y);
	  density = energy = -1e20;

#include "libinit.cx" // here we load the initial conditions

	  Tau = 2.0*M_PI/(v_azimuth+OMEGAFRAME);
	  rho0 = density; //initial condition
	  vrad0 = vrad;
	  vtheta0 = v_azimuth;

	  //lambda = 30*fabs(ramp*ramp)/Tau*dt;
	  lambda = 2200*fabs(ramp*ramp*ramp)/Tau*dt;

    //if(fw->Fluid->next == NULL)
    dens[m] = (dens[m]+lambda*rho0)/(1.+lambda);


	  vel[_RAD_][m] = (vel[_RAD_][m]+lambda*vrad0)/(1.+lambda);
	  vel[_AZIM_][m] = (vel[_AZIM_][m]+lambda*vtheta0)/(1.+lambda);
	  if (_COLAT_ < NDIM)
	    vel[_COLAT_][m] = vel[_COLAT_][m]/(1.+lambda);
	}
      }
    }
  }
}


void ApplyStockholmBoundaryConditions2 (dt)
real dt;
{
  FluidWork *fw;
  real *center[3];
  real *vel[3], *dens, ramp=0.0, lambda, R_in, R_out;
  real Tau=0.0, rho0, vrad0, vtheta0, Rmin, Rmax;
  real margin, highmargin, Colmin, Colmax, damp_zone_min, damp_zone_max;
  long i, j, k, m, gncell[3], stride[3], dim;
  boolean correct;
  real radius, azimuth, colatitude, x, y, z;
  real vx, vy, vz, vrad, v_azimuth, v_colatitude;
  real Bx, By, Bz, Brad, B_azimuth, B_colatitude;
  real density, energy, axial_radius, a2m, a2p;
  real interm1, interm2, interm3; /* Work variables for user convenience */
  margin = 0.01*PERCENTMARGIN;
  highmargin = 0.01*HIGHPERCENTMARGIN;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  SetFluidProperties (fw->Fluid);
  for (dim = 0; dim < 3; dim++) {
    center[dim] = fw->desc->Center[dim];
    vel[dim] = fw->Velocity[dim];
  }
  Rmin = corner_min0[_RAD_]; //inner domain boundary
  Rmax = corner_max0[_RAD_]; //inner domain boundary
  Colmin = corner_min0[_COLAT_]; //inner domain boundary
  Colmax = corner_max0[_COLAT_]; //inner domain boundary
  damp_zone_min = margin*fabs(Colmin-M_PI/2.0);
  if (ASPECTRATIO < damp_zone_min) damp_zone_min = ASPECTRATIO;
  damp_zone_max=margin*fabs(Colmax-M_PI/2.0);
  if (ASPECTRATIO < damp_zone_max) damp_zone_max = ASPECTRATIO;
  R_in = Rmin*(1.+margin); // boundary of the wave killing zone, at inner boundary
  R_out= Rmax/(1.+margin); // boundary of the wave killing zone, at outer boundary

  dens = fw->Density;

  for (k = 0; k < gncell[2]; k++) {
    for (j = 0; j < gncell[1]; j++) {
      for (i = 0; i < gncell[0]; i++) {
	       m = i*stride[0]+j*stride[1]+k*stride[2];
	       radius = center[_RAD_][m];
	       colatitude = center[_COLAT_][m];
	       correct = NO;
	       x = y = 0.0;
	       if (radius < R_in) { //we are in the inner damping zone
	          correct = YES; //here we need to damp
	           x = (R_in-radius)/(R_in-Rmin); //rampx
	       }

         if (radius > R_out) {//we are in the outer damping zone
	          correct = YES; //here we need to damp
	          x = (radius-R_out)/(Rmax-R_out);//rampx
	        }
          x *= x; //Parabolic ramp as in De Val Borro et al (2006)

         if (__SPHERICAL && (_COLAT_ < NDIM)) {
	          if (fabs(Colmin - colatitude) < damp_zone_min) {
	             correct = YES;
	             y = 1. - (colatitude - Colmin) / damp_zone_min;
	          }

	          if (fabs(Colmax - colatitude) < damp_zone_max) {
	             correct = YES;
	             y = 1. - (Colmax - colatitude) / damp_zone_max;
	          }
            y *= y;
	        }

          if (correct == TRUE) {
             ramp = x+y;
	           //ramp = 1.-(1.-x)*(1.-y);
	           density = energy = -1e20;

             #include "libinit.cx" // here we load the initial conditions

             Tau = 2.0*M_PI/(v_azimuth+OMEGAFRAME); //local orbital period
             Tau *= 1.0;
	           rho0 = density; //initial condition
	           vrad0 = vrad;
	           vtheta0 = v_azimuth;

             lambda = Tau/ramp;

             //if(fw->Fluid->next == NULL)
             dens[m] = (dens[m]*lambda+rho0*dt)/(dt+lambda);
             vel[_RAD_][m] = (vel[_RAD_][m]*lambda+vrad0*dt)/(dt+lambda);
             vel[_AZIM_][m] = (vel[_AZIM_][m]*lambda+vtheta0*dt)/(dt+lambda);
             if (_COLAT_ < NDIM)
         	      vel[_COLAT_][m] = vel[_COLAT_][m]*lambda/(dt+lambda);

	        }
      }
    }
  }
}
