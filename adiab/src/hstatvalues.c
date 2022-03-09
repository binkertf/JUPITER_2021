/** \file hstatvalues.c

    Functions used to integrate equilibrium quantities on interfaces,
    and azimuthal fluxes expressions.  Contains functions that return any
    equilibrium quantity at an arbitrary location, as well as functions
    that integrate a given field (ie, function) on an interface, in 2D or
    3D, in any geometry. Finally, contains the functions needed by the
    meridian flux corrections, which give the mass or momentum flux in the
    inertial or rotating frame, or in a frame that rotates solidly at rate
    OMEGAFRAME.

*/

#include "jupiter.h"
#include "init.h"

real EnergyTotFluxInertial (point)
     real point[3];
{
  real radius, azimuth, colatitude;
  real z, r;
  r = radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  if (__SPHERICAL)
    r = sin(colatitude) * radius;
  return IC_Energy_tot(point)*r*(IC_Vazimuth(point)+OMEGAFRAME);
}

real MassFluxInertial (point)
     real point[3];
{
  real radius, azimuth, colatitude;
  real z, r;
  r = radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  if (__SPHERICAL)
    r = sin(colatitude) * radius;
  return IC_Density(point)*r*(IC_Vazimuth(point)+OMEGAFRAME);
}

real MassFluxSolidRotation (point)
     real point[3];
{
  real radius, azimuth, colatitude;
  real z, r;
  r = radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  if (__SPHERICAL)
    r = sin(colatitude) * radius;
  return IC_Density(point)*r*OMEGAFRAME;
}

real EnergyTotFluxSolidRotation (point)
     real point[3];
{
  real radius, azimuth, colatitude;
  real z, r;
  r = radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  if (__SPHERICAL)
    r = sin(colatitude) * radius;
  return IC_Energy_tot(point)*r*OMEGAFRAME;
}

real MomentumFluxInertial (point)
     real point[3];
{
  real radius, azimuth, colatitude;
  real z, r;
  r = radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  if (__SPHERICAL)
    r = sin(colatitude) * radius;
  if (Isothermal)
    return IC_Density(point)*r*(r*r*(IC_Vazimuth(point)+OMEGAFRAME)*(IC_Vazimuth(point)+OMEGAFRAME) \
				+IC_Energy(point));
  else
    return IC_Density(point)*r*(r*r*(IC_Vazimuth(point)+OMEGAFRAME)*(IC_Vazimuth(point)+OMEGAFRAME) \
				+IC_Energy(point)*(GetGamma()-1.0)/IC_Density(point));
}

real MomentumFluxSolidRotation (point)
     real point[3];
{
  real radius, azimuth, colatitude;
  real z, r;
  r = radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  if (__SPHERICAL)
    r = sin(colatitude) * radius;
  return IC_Density(point)*r*(r*r*OMEGAFRAME*(IC_Vazimuth(point)+OMEGAFRAME));
}

real IC_Average_Simpson (begin_point, end_point, function, steps)
     real begin_point[3], end_point[3];
     real (*function)();
     long steps;
{
  real dx[3], point[3];
  long i, k;
  real sum=0.0, value;
  if (steps < 2) pError ("too small number of integration steps in IC_Average, hstatvalues.c");
  steps *= 2; 			/* We make sure 'steps' is even */
  for (k = 0; k < 3; k++) {
    dx[k] = (end_point[k]-begin_point[k])/(real)steps;
    point[k] = begin_point[k];
  }
  for (i = 0; i <= steps; i++) { /* Simpson's rule */
    value = function (point);
    if ((i == 0) || (i == steps))
      sum += value;
    else {
      if (i % 2 == 0) 
	sum += 2.*value;
      else
	sum += 4.*value;
    }
    for (k = 0; k < 3; k++) 
      point[k] += dx[k];
  }
  return sum/(3.0*(real)steps);
}

real IC_Average (begin_point, end_point, function, steps)
     real begin_point[3], end_point[3];
     real (*function)();
     long steps;
{
  real dx[3], point[3];
  long i, k;
  real sum=0.0, value;
  for (k = 0; k < 3; k++) {
    dx[k] = (end_point[k]-begin_point[k])/(real)steps;
    point[k] = begin_point[k];
  }
  for (i = 0; i <= steps; i++) { /* Trapezium rule */
    value = function (point);
    if ((i == 0) || (i == steps))
      sum += .5*value;
    else
      sum += value;
    for (k = 0; k < 3; k++)
      point[k] += dx[k];
  }
  return sum/((real)steps);
}

real IC_2D_Mean (begin_point, end_point, function, steps, dim)
     real begin_point[3], end_point[3];
     real (*function)();
     long steps, dim;		/* dim: dim perp to surface of integration */
{
  real point[3];
  long i, j, dim1, dim2;
  real sum=0.0, value, b, a;
  real am, ap, bm, bp, da, db;
  real amin, amax, bmin, bmax;
  real radius, colatitude;
  real area, total_area=0.0;
  radius = begin_point[0];
  colatitude = begin_point[2];
  dim1 = (dim == 0);
  dim2 = 2 - (dim == 2);
  cpTriplet (begin_point, point);
  amin = begin_point[dim1];
  amax = end_point[dim1];
  bmin = begin_point[dim2];
  bmax = end_point[dim2];
  da = (amax-amin)/(real)steps;
  db = (bmax-bmin)/(real)steps;
  for (i = 0; i <= steps; i++) { /* Trapezium rule */
    a = amin+(real)i*da;
    point[dim1] = a;
    am = a-.5*da;
    ap = a+.5*da;
    if (i == 0) am = amin;
    if (i == steps) ap = amax;
    for (j = 0; j <= steps; j++) {
      b = bmin+(real)j*db;
      point[dim2] = b;
      value = function (point);
      bm = b-.5*db;
      bp = b+.5*db;
      if (j == 0) bm = bmin;
      if (j == steps) bp = bmax;
      area = (ap-am)*(bp-bm);

      if (__CYLINDRICAL && (dim == 0)) /* radius in cylindrical coordinates */
	area *= radius;		

      if (__CYLINDRICAL && (dim == 2)) /* z in cylindrical coordinates */
	area = .5*(ap*ap-am*am)*(bp-bm);

      if (__SPHERICAL && (dim == 0)) /* radius in spherical coordinates */
	area = radius*radius*(ap-am)*(cos(bm)-cos(bp));

      if (__SPHERICAL && (dim == 1)) /* azimuth in spherical coordinates */
	area = (ap*ap-am*am)*.5*(bp-bm);

      if (__SPHERICAL && (dim == 2)) /* colatitude in spherical coordinates */
	area = .5*(ap*ap-am*am)*(bp-bm)*sin(colatitude);

      sum += area*value;
      total_area += area;
    }
  }
  return sum/total_area;
}

real IC_Corr (begin_point, end_point, function, steps, mid_point)
     real begin_point[3], end_point[3], mid_point[3];
     real (*function)();
     long steps;
{
  return IC_Average (begin_point, end_point, function, steps) / function(mid_point);
}

real IC_Pressure (point)
     real point[3];
{
  if(Isothermal)
    return IC_Density(point)*IC_Energy(point);
  else
    return IC_Energy(point)*(GetGamma()-1.);
}

real IC_PressureByRad (point)
     real point[3];
{
  real radius;
  radius = point[0];
  if(Isothermal)
    return IC_Density(point)*IC_Energy(point)*radius;
  else
    return IC_Energy(point)*(GetGamma()-1.)*radius;
}

real IC_Density (point)
     real point[3];
{
  return GetInitialValue (point, _density_);
}

real IC_Vrad (point)
     real point[3];
{
  return GetInitialValue (point, _vrad_);
}

real IC_Vazimuth (point)
     real point[3];
{
  return GetInitialValue (point, _vazimuth_);
}

real IC_Vcolatitude (point)
     real point[3];
{
  return GetInitialValue (point, _vcolatitude_);
}

real IC_Energy (point)
     real point[3];
{
  return GetInitialValue (point, _energy_);
}

real IC_Energy_tot (point)
     real point[3];
{
  return GetInitialValue (point, _tot_energy_);
}

real GetInitialValue (point, type)
     real point[3];
     long type;
{
  real radius, azimuth, colatitude;
  real x=0.0, y=0.0, z, value;
  real density, energy, vrad, v_azimuth, v_colatitude, vx, vy, vz;
  real interm1, interm2, interm3;
  real axial_radius;
  radius = point[0];
  azimuth = point[1];
  z = colatitude = point[2];
  density = energy = -1e20;
#include "libinit.cx"
  value = density;
  if (type == _vrad_)
    value = vrad;
  if (type == _vazimuth_)
    value = v_azimuth;
  if (type == _vcolatitude_)
    value = v_colatitude;
  if (type == _energy_)
    value = energy;
  if (type == _tot_energy_)
    value = energy+.5*pow(v_azimuth*radius,2.0)*density;
  return value;
}
