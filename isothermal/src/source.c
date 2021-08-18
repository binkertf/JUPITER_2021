#include "jupiter.h"

void Source (dt)
     real dt;
{
  FluidWork *fw;
  real *center[3];
  real *vel[3], source, *scor=NULL;
  long i, j, k, gncell[3], stride[3], dim;
  long strides[3];
  long m, ms=0, l;
  real metric;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  for (dim = 0; dim < 3; dim++) {
    center[dim] = fw->desc->Center[dim];
    vel[dim] = fw->Velocity[dim];
  }
  for (l = 0; l < NDIM; l++) {
    if (CorrHydroStat[l]) {
      scor = fw->Fluid->Source_corr[l]->field;
      for (i = 0; i < 3; i++)
	strides[i] = fw->Fluid->Rho_eq_c->stride[i];
    }
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
	for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  if (CorrHydroStat[l])
	    ms = i*strides[0]+j*strides[1]+k*strides[2];
	  metric = 1.0;
	  if (__CYLINDRICAL) {
	    if (l == _AZIM_) metric = center[_RAD_][m];
	  }
	  if (__SPHERICAL) {
	    if (l == _AZIM_) metric = center[_RAD_][m]*sin(center[_COLAT_][m]);
	    if (l == _COLAT_) metric = center[_RAD_][m];
	  }
	  source = fw->SourceVelocity[l][m];
	  if (CorrHydroStat[l])
	    source += scor[ms];
	  vel[l][m] += source*dt/metric;
	  /* The 'metric' correction is due to the fact that we are
	     possibly dealing with an angular velocity here. There is
	     no correction (metric=1.0) in the linear case*/

	}
      }
    }
  }
  if (fabs(DRIFTVELOCITY) > 1e-10)
    DiskRadialDrift(dt);
}


//------------------------------------------------------------------------------
static Beam beam, beam3,beam4;

void DustDiffusion(dt)
    real dt;
{
  long i, j, k, size[3],size2[3], dim, ip1, ip2; // variable declaration
  for (i = 0; i < 3; i++)
    size[i] = CurrentFluidPatch->desc->ncell[i];
  //printf("%s\n",SecondaryFluidPatch->Fluid->Name );
  //FillSources (PREDICT, EVERYWHERE);
  //FillSlopes ();
  for (dim = 0; dim < NDIM; dim++) { // loop over all dimensions (dim) -> RAD, COLAT, AZIM
    ip1 = (dim == 0);
    ip2 = 2-(dim == 2);
    AllocSecondaryBeam_PLM (&beam4, SecondaryFluidPatch->desc->gncell[dim]); //gas beam
    AllocTertiaryBeam_PLM (&beam3, CurrentFluidPatch->desc->gncell[dim]); // allocate beam for the dust fluid

    for (k = 0; k < size[ip2]; k++) {
      for (j = 0; j < size[ip1]; j++) {
        FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam3);
        FillBeam2 (dim, j+Nghost[ip1], k+Nghost[ip2], &beam4);
        gfo(&beam3, dt); // I am not sure if ths is needed here
        plm(&beam4, dt);

        Compute_Fluxes_Diffusion(&beam3,&beam4,dt); // compute the diffusion fluxes for that dim
        FillDiffFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam3); // fill fluxes for that dim
      }
    }

  }
    MassUpdate (dt); //dust diffusion as in Weber et al. (2019) without angular momentum conservation



}
