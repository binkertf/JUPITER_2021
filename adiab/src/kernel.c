  #include "jupiter.h"

static Beam beam, beam2;

void HydroKernel (dt)
     real dt;
{
  long i, j, k, size[3], dim, ip1, ip2, lev;
  lev = CurrentFluidPatch->desc->level;
  for (i = 0; i < 3; i++)
    size[i] = CurrentFluidPatch->desc->ncell[i];
  if (Stellar) //i.e. non isothermal
    ComputeQplus(); // viscous heating computation
  if (mGFO){ 
    JUP_SAFE(FillSources (PREDICT, EVERYWHERE));
  }
  if (mPLM){ 
    JUP_SAFE(FillSources (PREDICT, EVERYWHERE));
    JUP_SAFE(FillSlopes ());
  }
  if (mMUSCL) {
    JUP_SAFE(FillSources_Predict());
    if (!Isothermal){
      Predictor_adiab (dt);
    }else{
      Predictor (dt); //Not implemented correctly in isothermal...
      //prs_error ("ERROR: Isothermal + MUSCL not implemented correctly. Better use PLM or GFO method.");
    }
  }
  if (!Isothermal)
    FillEnergyTot();
  for (dim = 0; dim < NDIM; dim++) { /* For each dimension */
    ip1 = (dim == 0);
    ip2 = 2-(dim == 2);
    JUP_SAFE(AllocBeam_PLM (&beam, CurrentFluidPatch->desc->gncell[dim]));
    for (k = 0; k < size[ip2]; k++) {
      for (j = 0; j < size[ip1]; j++) {
	JUP_SAFE(FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	/* Scan a beam of the active mesh */
  JUP_SAFE(__Prepare_Riemann_States (&beam, dt));
	/* which is used to prepare the Riemann States */
	JUP_SAFE(__Compute_Fluxes (&beam, dt)); 
	/* The Riemann solver is then called and the fluxes evaluated */
	JUP_SAFE(FillFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	/* and fluxes are stored for that dim */
      }
    }
  }
  JUP_SAFE(ConservativeUpdate (dt));
  JUP_SAFE(DensFloor());
  /* and the conservative update is performed, together */
  /* with a face flux monitoring */
  JUP_SAFE(PressureCorrection (dt));
  /* Needs to be done *before* the source filling in SPHERICAL */
  JUP_SAFE(FillSources_geom_col (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  if (lev < HIGHRESLEVEL) {
    if (!Isothermal)
      JUP_SAFE(EnergyCorrection (dt));
  }
  /* Leave this after source step in order not to interfere with the
     divergence evaluation performed earlier. */
  JUP_SAFE(FillSources_geom_rad (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  JUP_SAFE(FillSources_pot (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  if (lev >= HIGHRESLEVEL) {
    if (!Isothermal)
      JUP_SAFE(EnergyCorrection2 (dt));
  }
  /* Apply source terms (potential gradient, centrifugal force) */

  if (KEPLERIAN && !NoStockholm)
    ApplyStockholmBoundaryConditionsDust (dt);
      //ApplyStockholmBoundaryConditions (dt);

    if (Stellar) {
    JUP_SAFE (RT_main(dt)); /* only apply to the gas fluid*/
    //if there is a dust fluid, fill dust temperature after gas temperature has been computed
    if (NbFluids>1){
      FillDust();
    }
  }
}

//######################################################################

//######################################################################


void DustKernel (dt)
      real dt;
{
  long i, j, k, size[3], dim, ip1, ip2, lev;
  lev = CurrentFluidPatch->desc->level;
  for (i = 0; i < 3; i++)
    size[i] = CurrentFluidPatch->desc->ncell[i];
  //JUP_SAFE(FillSources (PREDICT, EVERYWHERE));
  //JUP_SAFE(FillSources_Predict());
  //JUP_SAFE(FillSlopes ());
  //Predictor (dt); 
  for (dim = 0; dim < NDIM; dim++) { /* For each dimension */
    ip1 = (dim == 0);
    ip2 = 2-(dim == 2);
    JUP_SAFE(AllocBeam_PLM (&beam, CurrentFluidPatch->desc->gncell[dim]));
    for (k = 0; k < size[ip2]; k++) {
      for (j = 0; j < size[ip1]; j++) {
	      JUP_SAFE(FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	      /* Scan a beam of the active mesh */
        if (!Isothermal){
          JUP_SAFE(gfo_adiab(&beam, dt));
        }else{
          JUP_SAFE(gfo(&beam, dt));
        }
	      /* which is used to prepare the Riemann States */
	      JUP_SAFE(__Compute_Fluxes_pressureless(&beam, dt));
        /* The Riemann solver is then called and the fluxes evaluated */
	      JUP_SAFE(FillFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	      /* and fluxes are stored for that dim */
      }
    }
  }

  //if (DUSTDIFF == YES){
    //JUP_SAFE(DustDiffusion (dt));
  //}

  JUP_SAFE(ConservativeDustUpdate (dt));
  JUP_SAFE(DustDensFloor());

  /* and the conservative update is performed, together */
  /* with a face flux monitoring */
  //JUP_SAFE(PressureCorrection (dt));
  /* Needs to be done *before* the source filling in SPHERICAL */
  JUP_SAFE(FillSources_geom_col (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  //if (lev < HIGHRESLEVEL) {
    //if (!Isothermal)
      //JUP_SAFE(EnergyCorrection (dt));
  //}
  /* Leave this after source step in order not to interfere with the
     divergence evaluation performed earlier. */
  JUP_SAFE(FillSources_geom_rad (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  JUP_SAFE(FillSources_pot (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  //if (lev >= HIGHRESLEVEL) {
    //if (!Isothermal)
      //JUP_SAFE(EnergyCorrection2 (dt));
  //}
  /* Apply source terms (potential gradient, centrifugal force) */
  if (KEPLERIAN && !NoStockholm) {
    ApplyStockholmBoundaryConditionsDust (dt);
  }

}


//######################################################################

//######################################################################

void DustDiffPresKernel (dt)
      real dt;
{
  long i, j, k, size[3], dim, ip1, ip2, lev;
  lev = CurrentFluidPatch->desc->level;
  for (i = 0; i < 3; i++)
    size[i] = CurrentFluidPatch->desc->ncell[i];
  JUP_SAFE(FillSources (PREDICT, EVERYWHERE));
  //JUP_SAFE(FillSources_Predict());
  //JUP_SAFE(FillSlopes ());
  //Predictor (dt); 
  for (dim = 0; dim < NDIM; dim++) { /* For each dimension */
    ip1 = (dim == 0);
    ip2 = 2-(dim == 2);
    JUP_SAFE(AllocBeam_PLM (&beam, CurrentFluidPatch->desc->gncell[dim]));
    for (k = 0; k < size[ip2]; k++) {
      for (j = 0; j < size[ip1]; j++) {
	      JUP_SAFE(FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	      /* Scan a beam of the active mesh */
        JUP_SAFE(__Prepare_Riemann_States (&beam, dt));
	      /* which is used to prepare the Riemann States */
	      JUP_SAFE(Compute_Fluxes_Iso(&beam, dt));
        /* The Riemann solver is then called and the fluxes evaluated */
	      JUP_SAFE(FillFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	      /* and fluxes are stored for that dim */
      }
    }
  }

  JUP_SAFE(ConservativeUpdate (dt));
  JUP_SAFE(DustDensFloor());

  /* and the conservative update is performed, together */
  /* with a face flux monitoring */
  //JUP_SAFE(PressureCorrection (dt));
  /* Needs to be done *before* the source filling in SPHERICAL */
  JUP_SAFE(FillSources_geom_col (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  /* Leave this after source step in order not to interfere with the
     divergence evaluation performed earlier. */
  JUP_SAFE(FillSources_geom_rad (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  JUP_SAFE(FillSources_pot (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));

  JUP_SAFE(FillSources_diff (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));



  /* Apply source terms (potential gradient, centrifugal force) */
  if (KEPLERIAN && !NoStockholm) {
    ApplyStockholmBoundaryConditionsDust (dt);
  }

}


//######################################################################

//######################################################################




void DustDiffKernel (dt)
      real dt;
{
  long i, j, k, size[3], dim, ip1, ip2, lev;
  lev = CurrentFluidPatch->desc->level;
  for (i = 0; i < 3; i++)
    size[i] = CurrentFluidPatch->desc->ncell[i];
  JUP_SAFE(FillSources (PREDICT, EVERYWHERE));
  //JUP_SAFE(FillSlopes ());
  //mMUSCL = NO;
  //JUP_SAFE(FillSources_Predict());
  //mMUSCL = YES;

  for (dim = 0; dim < NDIM; dim++) { /* For each dimension */
    ip1 = (dim == 0);
    ip2 = 2-(dim == 2);
    JUP_SAFE(AllocBeam_PLM (&beam, CurrentFluidPatch->desc->gncell[dim]));
    JUP_SAFE(AllocSecondaryBeam_PLM (&beam2, CurrentFluidPatch->desc->gncell[dim]));
    for (k = 0; k < size[ip2]; k++) {
      for (j = 0; j < size[ip1]; j++) {
	JUP_SAFE(FillBeam (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
  JUP_SAFE(FillBeam2 (dim, j+Nghost[ip1], k+Nghost[ip2], &beam2));
	/* Scan a beam of the active mesh */
	if (!Isothermal){
     JUP_SAFE(gfo_adiab(&beam, dt));
  }else{
    JUP_SAFE(gfo(&beam, dt));
  }
	/* which is used to prepare the Riemann States */
	JUP_SAFE(__Compute_Fluxes_pressureless(&beam, dt));
  JUP_SAFE(Compute_Fluxes_Diffusion(&beam, &beam2, dt));
	/* The Riemann solver is then called and the fluxes evaluated */
	JUP_SAFE(FillFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
  JUP_SAFE(FillDiffFluxes (dim, j+Nghost[ip1], k+Nghost[ip2], &beam));
	/* and fluxes are stored for that dim */
      }
    }
  }


  DiffusionUpdate (dt);
  JUP_SAFE(ConservativeDustUpdate (dt));
  JUP_SAFE(DustDensFloor());



  /* and the conservative update is performed, together */
  /* with a face flux monitoring */
  //JUP_SAFE(PressureCorrection (dt));
  /* Needs to be done *before* the source filling in SPHERICAL */
  JUP_SAFE(FillSources_geom_col (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  //if (lev < HIGHRESLEVEL) {
    //if (!Isothermal)
      //JUP_SAFE(EnergyCorrection (dt));
  //}
  /* Leave this after source step in order not to interfere with the
     divergence evaluation performed earlier. */
  JUP_SAFE(FillSources_geom_rad (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  JUP_SAFE(FillSources_pot (UPDATE, EVERYWHERE));
  JUP_SAFE(Source (dt));
  //if (lev >= HIGHRESLEVEL) {
    //if (!Isothermal)
      //JUP_SAFE(EnergyCorrection2 (dt));
  //}
  /* Apply source terms (potential gradient, centrifugal force) */





  if (KEPLERIAN && !NoStockholm) {
    ApplyStockholmBoundaryConditionsDust (dt);
  }

}
