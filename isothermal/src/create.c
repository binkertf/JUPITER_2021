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
  real *MomCorrRatioL[3][3], *MomCorrRatioR[3][3];
  real *Metric[3][2], *InvMetric[3][2];
  real frac;
  long size, MeridianSize, count=0;
  FluidPatch *previousFluid = NULL;
  FluidPatch *fluid;
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
	  Center[CoordNb[0]][l] = .5*(rmin+rmax);
	  Center[CoordNb[1]][l] = .5*(phimin+phimax);
	  Center[CoordNb[2]][l] = .5*(thetamin+thetamax);	/* approximate only */
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
				/* Z coordinate */
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
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      MomentumCorrection[i][j] = FALSE;
      MomCorrRatioL[i][j] = NULL;
      MomCorrRatioR[i][j] = NULL;
    }
  }
  if (__CYLINDRICAL) {
    MomentumCorrection[_RAD_][_AZIM_] = TRUE;
    MomCorrRatioL[_RAD_][_AZIM_] = prs_malloc (gncell[_RAD_]*sizeof(real));
    MomCorrRatioR[_RAD_][_AZIM_] = prs_malloc (gncell[_RAD_]*sizeof(real));
    for (i = 0; i < gncell[_RAD_]; i++) {
      l = i * stride[_RAD_];
      MomCorrRatioL[_RAD_][_AZIM_][i] = Center[_RAD_][l]/Edges[_RAD_][i+1];
      MomCorrRatioR[_RAD_][_AZIM_][i] = Center[_RAD_][l]/Edges[_RAD_][i];
    }
  }
  if (__SPHERICAL) {
    MomentumCorrection[_RAD_][_AZIM_] = TRUE;
    MomCorrRatioL[_RAD_][_AZIM_] = prs_malloc (gncell[_RAD_]*sizeof(real));
    MomCorrRatioR[_RAD_][_AZIM_] = prs_malloc (gncell[_RAD_]*sizeof(real));
    for (i = 0; i < gncell[_RAD_]; i++) {
      l = i * stride[_RAD_];
      MomCorrRatioL[_RAD_][_AZIM_][i] = Center[_RAD_][l]/Edges[_RAD_][i+1];
      MomCorrRatioR[_RAD_][_AZIM_][i] = Center[_RAD_][l]/Edges[_RAD_][i];
    }
    MomentumCorrection[_COLAT_][_AZIM_] = TRUE;
    MomCorrRatioL[_COLAT_][_AZIM_] = prs_malloc (gncell[_COLAT_]*sizeof(real));
    MomCorrRatioR[_COLAT_][_AZIM_] = prs_malloc (gncell[_COLAT_]*sizeof(real));
    for (i = 0; i < gncell[_COLAT_]; i++) {
      l = i * stride[_COLAT_];
      MomCorrRatioL[_COLAT_][_AZIM_][i] = sin(Center[_COLAT_][l])/sin(Edges[_COLAT_][i+1]);
      MomCorrRatioR[_COLAT_][_AZIM_][i] = sin(Center[_COLAT_][l])/sin(Edges[_COLAT_][i]);
    }
    MomentumCorrection[_RAD_][_COLAT_] = TRUE;
    MomCorrRatioL[_RAD_][_COLAT_] = prs_malloc (gncell[_RAD_]*sizeof(real));
    MomCorrRatioR[_RAD_][_COLAT_] = prs_malloc (gncell[_RAD_]*sizeof(real));
    for (i = 0; i < gncell[_RAD_]; i++) {
      l = i * stride[_RAD_];
      MomCorrRatioL[_RAD_][_COLAT_][i] = Center[_RAD_][l]/Edges[_RAD_][i+1];
      MomCorrRatioR[_RAD_][_COLAT_][i] = Center[_RAD_][l]/Edges[_RAD_][i];
    }
  }
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      desc->MomCorrRatioL[i][j] = MomCorrRatioL[i][j];
      desc->MomCorrRatioR[i][j] = MomCorrRatioR[i][j];
    }
  }
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
