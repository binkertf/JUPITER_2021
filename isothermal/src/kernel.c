#include "jupiter.h"

static Beam beam;

void HydroKernel (dt)
     real dt;
{
  long i, j, k, size[3],size2[3], dim, ip1, ip2;
  for (i = 0; i < 3; i++)
    size[i] = CurrentFluidPatch->desc->ncell[i];

  FillSources (PREDICT, EVERYWHERE);
  FillSlopes ();
  if (mMUSCL) Predictor (dt);
  for (dim = 0; dim < NDIM; dim++) { /* For each dimension */
    ip1 = (dim == 0);
    ip2 = 2-(dim == 2);
    AllocBeam_PLM (&beam, CurrentFluidPatch->desc->gncell[dim]); //

    for (k = 0; k < size[ip2]; k++) {
      for (j = 0; j < size[ip1]; j++) {
				/* Scan a beam of the active mesh */
         if (CurrentFluidPatch->Fluid->next == NULL){ //gas fluid
           FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam);
            /* gas fluid computed here*/
            __Prepare_Riemann_States (&beam, dt);
				    /* is used to predict the Riemann States at the interface (with PLM, GFO or MUSCL method) */
	          __Compute_Fluxes (&beam, dt);
				     /* The Riemann solver is then called and the fluxes evaluated */
         }else{ //dust fluid
           FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam);
           /* dust fluid is computed here*/
            gfo(&beam, dt);
           /* the dust fluxes are computed here*/
            __Compute_Fluxes_pressureless (&beam, dt);
         }
	FillFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam);
				/* and fluxes are stored for that dim */
        /* the linear fluxes are also substituted for angular fluxes if applicable (cylindrical or spherical)*/
      }
    }
  }




  // Diffusion module
  if ((CurrentFluidPatch->Fluid->next != NULL)&&(DUSTDIFF == YES)){//we apply the diffusion module to the dust fluid
    DustDiffusion(dt);
  }


/*  if (VISCOSITY > 1e-16)
    ApplyViscousStress (dt);*/ // this should go in flux.c before the momenta update...
  CorrectFluxesFromFinerLevel ();
				/* the fluxes are then corrected from finer level */
  ConservativeUpdate (dt);
				/* and the conservative update is performed, together */
				/* with a face flux monitoring */
  PressureCorrection (dt);
				/* Needs to be done *before* the source filling in SPHERICAL */
        //  FillSources (UPDATE, EVERYWHERE);

  /* Apply source terms (potential gradient, centrifugal force) */
  Source (dt);

  /* add density floor here */

real *e, *rho;
long m, sized;
sized  = CurrentFluidPatch->desc->gncell[0];
sized *= CurrentFluidPatch->desc->gncell[1];
sized *= CurrentFluidPatch->desc->gncell[2];
rho = CurrentFluidPatch->Density;
for (m=0; m < sized; m++) {
  if (rho[m] < DUSTDENSFLOOR) {
    rho[m] = DUSTDENSFLOOR;
  }
}


  if ((KEPLERIAN && !NoStockholm ) && (CurrentFluidPatch->Fluid->next == NULL)){
    ApplyStockholmBoundaryConditions2 (dt); //gas
  }
  if ((KEPLERIAN && !NoStockholm ) && (CurrentFluidPatch->Fluid->next != NULL)){
    ApplyStockholmBoundaryConditions2 (dt); //dust
  }


}
