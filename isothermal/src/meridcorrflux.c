/** \file meridcorrflux.c

  Evaluate the flux correction in azimuth needed for a wakeless nested
grid system. The function InitFluxCorrectionMeridian() evaluates the
coefficients by which the mass and momenta meridional fluxes (i.e., in
azimuth) must be multiplicated in order to yield resolution
independent fluxes (that is, the sum of two neighboring fluxes on a
fine mesh must coincide with the corresponding unique flux on the
coarse mesh). The fluxes need be evaluated in the inertial frame (so
as to avoid accuracy issues near corotation). The additivity
requirement above is enforced by performing an integral on the
interface with a number of steps that depend on the resolution
(namely, 2^(highest level+1-current level)), where "highest level" is
the highest refinement level met on the mesh in the ring whose
intersection with the coarsest mesh is the underlying current coarsest
zone. It is not the highest level over the whole mesh, as that could
be computationnally very expensive.

 */

#include "jupiter.h"

void InitFluxCorrectionMeridian (f)
     FluidPatch *f;
{
  long nrad, ncol, i, j, ir, ic, gncell[3], stride[3], *ml;
  real *enerc1, *enerc2, *massc1, *massc2, *momc1, *momc2;
  real corr_massflux1=0.0, corr_massflux2=0.0;
  real corr_momentumflux1=0.0, corr_momentumflux2=0.0;
  real corr_energyflux1=0.0, corr_energyflux2=0.0;
  /* pc : point at corner, pmid : point in the middle */
  real pcmin[3], pcmaxr[3], pcmaxc[3], pcmaxt[3], pmidr[3], pmidc[3], pmid[3];
  long nsteps, n1, n2, m;
  fflush (stdout);
  getgridsize (f->desc, gncell, stride);
  SetFluidProperties (f);
  nrad = gncell[_RAD_];
  ncol = gncell[_COLAT_];
  massc1 = f->MassFluxCorrection1;
  massc2 = f->MassFluxCorrection2;
  momc1 = f->MomentumFluxCorrection1;
  momc2 = f->MomentumFluxCorrection2;
  enerc1 = f->EnergyFluxCorrection1;
  enerc2 = f->EnergyFluxCorrection2;
  ml = f->desc->MaxLevMeridianProj;
  n1 = gncell[_RAD_];
  n2 = gncell[_COLAT_];
  nsteps = 2L<<(LevMax-f->desc->level);
  if (_COLAT_ < _RAD_) swapl (&n1, &n2);
  for (i = 0; i < n1; i++) {
    for (j = 0; j < n2; j++) {
      ir = (_RAD_ < _COLAT_ ? i : j); 
      ic = (_RAD_ < _COLAT_ ? j : i); 
      m = ir*stride[_RAD_]+ic*stride[_COLAT_];
      nsteps = 2L<<(1+ml[i+j*n1]-f->desc->level);
      pcmin[0] = f->desc->Edges[_RAD_][ir];
      pcmin[1] = pmid[1]= 0.0;
      pcmin[2] = f->desc->Edges[_COLAT_][ic];
      cpTriplet (pcmin, pcmaxr);
      cpTriplet (pcmin, pcmaxc);
      cpTriplet (pcmin, pmidr);
      cpTriplet (pcmin, pmidc);
      pcmaxr[0] = f->desc->Edges[_RAD_][ir+1];
      cpTriplet (pcmaxr, pcmaxt);
      pcmaxt[2] = f->desc->Edges[_COLAT_][ic+1];
      pcmaxc[2] = f->desc->Edges[_COLAT_][ic+1];
      pmid[0] = pmidr[0] = f->desc->Center[_RAD_][m];
      pmid[2] = pmidc[2] = f->desc->Center[_COLAT_][m];
      if (NDIM == 2) {
	corr_massflux1 = IC_Corr (pcmin,pcmaxr,MassFluxInertial,nsteps,pmidr);
	corr_massflux2 = IC_Corr (pcmin,pcmaxr,MassFluxSolidRotation,nsteps,pmidr);
	corr_momentumflux1 = IC_Corr (pcmin,pcmaxr,MomentumFluxInertial,nsteps,pmidr);
	corr_momentumflux2 = IC_Corr (pcmin,pcmaxr,MomentumFluxSolidRotation,nsteps,pmidr);
	if (!Isothermal) {	/* not correctly implemented yet cuidadin */
	  corr_energyflux1 = IC_Corr (pcmin,pcmaxr,MassFluxRotating,nsteps,pmidr);
	  corr_energyflux2 = IC_Corr (pcmin,pcmaxr,MassFluxRotating,nsteps,pmidr);
	}
      } else if (NDIM == 3) {
	corr_massflux1 = IC_2D_Mean (pcmin, pcmaxt, MassFluxInertial, nsteps, 1); /* 1: azimuthal */
	corr_massflux1 /= MassFluxInertial (pmid);
	corr_massflux2 = IC_2D_Mean (pcmin, pcmaxt, MassFluxSolidRotation, nsteps, 1); /* 1: azimuthal */
	corr_massflux2 /= MassFluxSolidRotation (pmid);
	corr_momentumflux1 = IC_2D_Mean (pcmin, pcmaxt, MomentumFluxInertial, nsteps, 1);
	corr_momentumflux1 /= MomentumFluxInertial (pmid);
	corr_momentumflux2 = IC_2D_Mean (pcmin, pcmaxt, MomentumFluxSolidRotation, nsteps, 1);
	corr_momentumflux2 /= MomentumFluxSolidRotation (pmid);
      }
      massc1[i+j*n1] = corr_massflux1;
      massc2[i+j*n1] = corr_massflux2;
      momc1[i+j*n1] = corr_momentumflux1;
      momc2[i+j*n1] = corr_momentumflux2;
      if (!Isothermal) {
	enerc1[i+j*n1] = corr_energyflux1;
	enerc2[i+j*n1] = corr_energyflux2;
      }
    }
  }
}
