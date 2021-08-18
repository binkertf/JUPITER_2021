/** \file hstatprep.c

  Allocate and affect the arrays used to store the equilibrium configuration.
  Prepare the arrays (possibly squeezed if the equilibrium state
  possesses some symmetry properties) that are used to store the
  information relative to the equilibrium state, and affect these
  arrays. The last function in this file is called by GlobalInit().
  It temporarily overwrites the InitCode and PotentialCode variables
  so as to fill the HD arrays with the equilibrium state.

*/

#include "jupiter.h"

extern char *CoordNames[];

SqueezedField *AllocateSqueezedArray (patch, name)
FluidPatch *patch;
char *name;
{
  long gncell[3], stride[3], i;
  long Size;
  SqueezedField *sf;
  getgridsize (patch->desc, gncell, stride);
  sf = (SqueezedField *)prs_malloc(sizeof(SqueezedField));
  sf->Name = (char *)prs_malloc(sizeof(char)*MAXLINELENGTH);
  strcpy (sf->Name, name);
  sf->desc = patch->desc;
  for (i = 0; i < 3; i++) {
    if (SqueezeDim[i])
      gncell[i] = 1;
    stride[i] = (i == 0 ? 1 : stride[i-1]*gncell[i-1]);
    sf->gncell[i] = gncell[i];
  }
  for (i = 0; i < 3; i++) {
    if (SqueezeDim[i])
      stride[i] = 0;
    sf->stride[i] = stride[i];
  }
  Size = gncell[0]*gncell[1]*gncell[2];
  sf->field = (real *)prs_malloc(sizeof(real)*Size);
  return sf;
}

void CreateSqueezedArrays (patch)
     FluidPatch *patch;
{
  long i;
  char coordname[MAXLINELENGTH], name[MAXLINELENGTH];
  patch->Rho_eq_c = AllocateSqueezedArray (patch, "HydroStatDensC");
  if (!Isothermal)
    patch->Ene_eq_c = AllocateSqueezedArray (patch, "HydroStatEnergyC");
  for (i = 0; i < 3; i++) {
    if (CorrHydroStat[i]) {
      strcpy (coordname, CoordNames[CoordType*3+InvCoordNb[i]]);

      sprintf (name, "HydroStatDensI_%s", coordname);
      patch->Rho_eq_i[i] = AllocateSqueezedArray (patch, name);

      if (!Isothermal) {
	sprintf (name, "HydroStatEnergyI_%s", coordname);
	patch->Ene_eq_i[i] = AllocateSqueezedArray (patch, name);
      }

      sprintf (name, "HydroCs2I_%s", coordname);
      patch->Cs2_i[i] = AllocateSqueezedArray (patch, name);

      sprintf (name, "HydroStatPresI_%s", coordname);
      patch->Pres_eq_i[i] = AllocateSqueezedArray (patch, name);

      sprintf (name, "HydroStatCorr_%s", coordname);
      patch->Source_corr[i] = AllocateSqueezedArray (patch, name);
    }
  }
}

void HydroStatPrepare ()
{
  long TempPotCode;
  long hydrostatnumb=99999L;
  boolean TempRestart;
  real TempGlobalDate;
  TempRestart = Restart;
  TempGlobalDate = GlobalDate;
  Restart = NO;
  GlobalDate = 0.0;
  //  TempInitCode = InitCode;
  InitMode = EQUILIBRIUM;
  TempPotCode = PotentialCode;
  //  findcodeinit (INITHYDROSTAT, &InitCode_Eq);
  findcodepot (POTENTIALHYDROSTAT, &PotCode_Eq);
  //InitCode = InitCode_Eq;
  PotentialCode = PotCode_Eq;
  InitWholeHierarchy (0L);
  /* The meshes now all contain what we consider to be the hydrostatic
     equilibrium */
  pInfo ("We output for further reference the equilibrium state\n");
  pInfo ("with output number %ld\n", hydrostatnumb);
  Write (hydrostatnumb);		/* We output it */
  /* We now create and initialize the squeezed arrays */
  ForAllPatches (CreateSqueezedArrays);
  ForAllPatches (SetCenteredDensity);
  ResetChronometer (3);
  ForAllPatches (SetInterfaceQuantities);
  ReadChronometer (3, "Interface quantities construction");
  /* and we initialize the source term correction */
  ForAllPatches (CorrectSourceTerm);
  WriteHydroStat (hydrostatnumb);
  /* Finally we have to set all HD arrays to zero before returning,
     since the Initialization function is incremental */
  ForAllPatches (ResetPatch);
  //InitCode = TempInitCode;
  InitMode = STANDARD;
  PotentialCode = TempPotCode;
  Restart = TempRestart;
  GlobalDate = TempGlobalDate;
}
