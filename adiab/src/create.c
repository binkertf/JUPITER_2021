#include "jupiter.h"

void FillCPUGrid (desc)
     tGrid_CPU *desc;
{
  long i, j, k, l, ip[2];
  long idx[3];
  long gncell[3], stride[3];
  real xmin, xmax, ymin, ymax, zmin, zmax;
  real rmin, rmax, phimin, phimax;
  real thetamin, thetamax, dl, coord;
  real volume=1.0, *InvVolume, **Center, **InterSurface, *Edges[3];
  real *Metric[3][2], *InvMetric[3][2];
  real frac;
  long size, MeridianSize, count=0;
  FluidPatch *previousFluid = NULL;
  FluidPatch *fluid;
  real *optical_depth;
  getgridsize (desc, gncell, stride);
  size = gncell[0]*gncell[1]*gncell[2];
  InvVolume = prs_malloc (size*sizeof(real));
  desc->InvVolume = InvVolume;
  InterSurface = multiple_alloc_1D (3, size);
  Center       = multiple_alloc_1D (3, size);
  for (i = 0; i < 3; i++) {	/* One needs to leave 3 here, instead
				   of NDIM, in case one has to refer
				   to a non-indifferent extra
				   dimension if NDIM<3 */
    ip[0] = (i == 0);
    ip[1] = 2 - (i == 2);
    desc->Center[i] = Center[i];
    desc->InterSurface[i] = InterSurface[i];
    Edges[i] = prs_malloc ((gncell[i]+1)*sizeof(real));
    for (l = 0; l < 2; l++) {
      Metric[i][l] = prs_malloc (gncell[ip[l]]*sizeof(real));
      desc->Metric[i][l] = Metric[i][l];
      InvMetric[i][l] = prs_malloc (gncell[ip[l]]*sizeof(real));
      desc->InvMetric[i][l] = InvMetric[i][l];
    }
    desc->Edges[i] = Edges[i];
  }
  for (j = 0; j < 3; j++) {
    for (i = 0; i <= gncell[j]; i++) {
      frac = (real)(i*desc->Parent->dn[j]+desc->gncorner_min[j])/\
	(real)(ncorner_max0[j]-ncorner_min0[j]);
      Edges[j][i] = frac*(corner_max0[j]-corner_min0[j])+corner_min0[j];
    }
  }
  optical_depth = prs_malloc (gncell[0]*gncell[2]*sizeof(real));
  desc->optical_depth = optical_depth;
  for (k = 0; k < gncell[2]; k++) {
    for (j = 0; j < gncell[1]; j++) {
      for (i = 0; i < gncell[0]; i++) {
	l = i*stride[0]+j*stride[1]+k*stride[2];
	idx[0] = i;
	idx[1] = j;
	idx[2] = k;
	switch (CoordType) {
	case CARTESIAN: 
	  xmin = Edges[CoordNb[0]][idx[CoordNb[0]]];
	  xmax = Edges[CoordNb[0]][idx[CoordNb[0]]+1];
	  ymin = Edges[CoordNb[1]][idx[CoordNb[1]]];
	  ymax = Edges[CoordNb[1]][idx[CoordNb[1]]+1];
	  zmin = Edges[CoordNb[2]][idx[CoordNb[2]]];
	  zmax = Edges[CoordNb[2]][idx[CoordNb[2]]+1];
	  Center[CoordNb[0]][l] = .5*(xmin+xmax);
	  Center[CoordNb[1]][l] = .5*(ymin+ymax);
	  Center[CoordNb[2]][l] = .5*(zmin+zmax);
	  InterSurface[CoordNb[0]][l] = InterSurface[CoordNb[1]][l] =\
	    InterSurface[CoordNb[2]][l] = volume = 1.0;
	  /* X coordinate */
	  dl = (_X_ < NDIM ? xmax-xmin : 1.0);
	  volume *= dl;
	  InterSurface[CoordNb[1]][l] *= dl;
	  InterSurface[CoordNb[2]][l] *= dl;
	  /* Y coordinate */
	  dl = (_Y_ < NDIM ? ymax-ymin : 1.0);
	  volume *= dl;
	  InterSurface[CoordNb[0]][l] *= dl;
	  InterSurface[CoordNb[2]][l] *= dl;
	  /* Z coordinate */
	  dl = (_Z_ < NDIM ? zmax-zmin : 1.0);
	  volume *= dl;
	  InterSurface[CoordNb[0]][l] *= dl;
	  InterSurface[CoordNb[1]][l] *= dl;
	  break;
	case CYLINDRICAL:
	  rmin =   Edges[CoordNb[0]][idx[CoordNb[0]]];
	  rmax =   Edges[CoordNb[0]][idx[CoordNb[0]]+1];
	  phimin = Edges[CoordNb[1]][idx[CoordNb[1]]];
	  phimax = Edges[CoordNb[1]][idx[CoordNb[1]]+1];
	  zmin =   Edges[CoordNb[2]][idx[CoordNb[2]]];
	  zmax =   Edges[CoordNb[2]][idx[CoordNb[2]]+1];
	  Center[CoordNb[0]][l] = .5*(rmin+rmax);
	  Center[CoordNb[1]][l] = .5*(phimin+phimax);
	  Center[CoordNb[2]][l] = .5*(zmin+zmax);
	  InterSurface[CoordNb[0]][l] = InterSurface[CoordNb[1]][l] =\
	    InterSurface[CoordNb[2]][l] = volume = 1.0;
	  /* Radial coordinate */
	  dl = (_RAD_ < NDIM ? .5*(rmax*rmax-rmin*rmin) : rmin);
	  volume *= dl;
	  InterSurface[CoordNb[1]][l] *= (_RAD_ < NDIM ? rmax-rmin : 1.0);
	  InterSurface[CoordNb[2]][l] *= dl;
	  /* Azimuth coordinate */
	  dl = (_AZIM_ < NDIM ? phimax-phimin : 1.0);
	  volume *= dl;
	  InterSurface[CoordNb[0]][l] *= dl * rmin;
	  InterSurface[CoordNb[2]][l] *= dl;
	  /* Z coordinate */
	  dl = (_Z_ < NDIM ? zmax-zmin : 1.0);
	  volume *= dl;
	  InterSurface[CoordNb[0]][l] *= dl;
	  InterSurface[CoordNb[1]][l] *= dl;
	  break;
	case SPHERICAL:
	  rmin     = Edges[CoordNb[0]][idx[CoordNb[0]]];
	  rmax     = Edges[CoordNb[0]][idx[CoordNb[0]]+1];
	  phimin   = Edges[CoordNb[1]][idx[CoordNb[1]]];
	  phimax   = Edges[CoordNb[1]][idx[CoordNb[1]]+1];
	  thetamin = Edges[CoordNb[2]][idx[CoordNb[2]]];
	  thetamax = Edges[CoordNb[2]][idx[CoordNb[2]]+1];
	  /* Intersurface is strictly homogoneous to length^2 */
	  Center[CoordNb[0]][l] = (2./3.)*(rmax*rmax*rmax-rmin*rmin*rmin)/(rmax*rmax-rmin*rmin);
	  Center[CoordNb[1]][l] = .5*(phimin+phimax);
	  if (.5*(thetamin+thetamax) < PI/2) 
	    Center[CoordNb[2]][l] = asin((cos(thetamin)-cos(thetamax))/(thetamax-thetamin));
	  else
	    Center[CoordNb[2]][l] = PI-asin((cos(thetamin)-cos(thetamax))/(thetamax-thetamin));
	  InterSurface[CoordNb[0]][l] = InterSurface[CoordNb[1]][l] =\
	    InterSurface[CoordNb[2]][l] = volume = 1.0;
	  /* Radial coordinate */
	  dl = (_RAD_ < NDIM ? ONETHIRD*(rmax*rmax*rmax-rmin*rmin*rmin) : rmin*rmin);
	  volume *= dl;
	  dl = (_RAD_ < NDIM ? .5*(rmax*rmax-rmin*rmin) : rmin);
	  InterSurface[CoordNb[1]][l] *= dl;
	  InterSurface[CoordNb[2]][l] *= dl;
	  /* Azimuth coordinate */
	  dl = (_AZIM_ < NDIM ? phimax-phimin : 1.0);
	  volume *= dl;
	  InterSurface[CoordNb[0]][l] *= dl * rmin * rmin;
	  InterSurface[CoordNb[2]][l] *= dl * sin(thetamin);
	  /* colatitude coordinate */
	  dl = (_COLAT_ < NDIM ? (cos(thetamin)-cos(thetamax)) : sin(thetamin));
	  volume *= dl;
	  InterSurface[CoordNb[0]][l] *= dl;
	  InterSurface[CoordNb[1]][l] *= (_COLAT_ < NDIM ? thetamax-thetamin : 1.0);
	  break;
	}
	InvVolume[l] = 1.0/volume;
      }
    }
  }
  for (i = 0; i < 3; i++) {	/* One needs to leave 3 here, instead
				   of NDIM, in case one has to refer
				   to a non-indifferent extra
				   dimension if NDIM<3 */
    ip[0] = (i == 0);
    ip[1] = 2 - (i == 2);
    for (l = 0; l < 2; l++) {
      for (j = 0; j < gncell[ip[l]]; j++) {
	coord = Center[ip[l]][j*stride[ip[l]]];
	if (__CARTESIAN)
	  coord = 1.0;
	if (__CYLINDRICAL && (!((i == _AZIM_) && (ip[l] == _RAD_))))
	  coord = 1.0;
	if (__SPHERICAL && ((i == _RAD_) || ((i == _COLAT_) && (ip[l] == _AZIM_))))
	  coord = 1.0;
	if (__SPHERICAL && (i == _AZIM_) && (ip[l] == _COLAT_))
	  coord = sin(coord);
	Metric[i][l][j] = coord;
	InvMetric[i][l][j] = 1./coord;
      }
    }
  }

  /* A comment on the Metric[3][2] array that should have been written many years ago... */
  /* In cartesian coordinates this matrix has just ones everywhere.
     In cylindrical coordinates (with COORDPERMUT set to 213) it has the following shape:

    [[r 1]
     [1 1]
     [1 1]]

 and in spherical coordinates (COORDPERMUT set to 213 also):

    [[r sin\theta]
     [1   1      ]
     [1   r      ]] */

  /* This metric matrix gives the length of an arc of a given
coordinate by simply multiplying, on a given line, the first element
by the second. Example: the length of an azimuthal arc of size d\phi
in cylindrical coordinates (azimuth == 1st coord ==> first line of
matrix) is given by: dl = d\phi * r * 1.

Another example: the length of a colatitude arc of size d\theta in
spherical coordinates is given by (3rd line, last coordinate...): 
dl = d\theta * 1 * r.

Another example: the length of an azimuthal arc of size d\phi in
spherical coordinates is given by (1st line, first coordinate...): 
dl = d\phi * r * sin(\theta).

The matrix is correctly evaluated regardless of the value of COORDPERMUT.
It should always feature in pair, involving the multiplication of the first
element by the second element. So we should always see a [0] ... [1] ...
in lines involving the metric.

The reason why we have to multiply the elements (instead of storing a single
array that is the result of the multiplication) is that we only need, this
way, to store 1D array, with a much lighter memory imprint.

The array InvMetric has just coefficients which are 1 by 1 the inverse of those
of Metric.
  */


  if (KEPLERIAN) {
    MeridianSize = gncell[_RAD_]*gncell[_COLAT_];
    desc->MaxLevMeridianProj      = prs_malloc (sizeof(long)*MeridianSize);
  }
  for (i = 0; i < NbFluids; i++) {
    desc->Fluid = CreateFluidPatch (desc, FluidName[i], InitCodeNames[i], InitCodeNamesEq[i]);
    desc->Fluid->next = previousFluid;
    if (previousFluid != NULL)
      previousFluid->prev = desc->Fluid;
    previousFluid = desc->Fluid;
  }
  fluid = desc->Fluid;
  count = 0;
  while (fluid != NULL) {
    fluid->FluidRank = count;
    count++;
    fluid = fluid->next;
  }
}
