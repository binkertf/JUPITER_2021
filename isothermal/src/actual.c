#include "jupiter.h"

static long MaxWorkSize=0, MaxWorkSize2=0;
static long AllocPhase, AllocPhase2;
static boolean NeverAllocated = TRUE,NeverAllocated2 = TRUE;
static long Length = 0, Length2 = 0;
static long NField = 0, NField2 = 0, FieldIndex = 0,FieldIndex2 = 0;
static real *FieldStart,  *FieldStart2;
static long PreviousLength=0, PreviousLength2=0;

void ReserveForWorkPatch (ptr)
     real **ptr;
{
  if (!AllocPhase)
    NField++;
  else {
    *ptr = FieldStart+Length*FieldIndex;
    FieldIndex++;
  }
}

void MakeCurrentFluidPatch ()
{
  long i, j;
  ReserveForWorkPatch (&(CurrentFluidPatch->Density));
  ReserveForWorkPatch (&(CurrentFluidPatch->Energy));
  ReserveForWorkPatch (&(CurrentFluidPatch->Potential));
  ReserveForWorkPatch (&(CurrentFluidPatch->SourceRhoPred));
  ReserveForWorkPatch (&(CurrentFluidPatch->Density_Pred));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Velocity[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Velocity_Pred[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->SourceVelocity[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Accel[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->SlopeDensity[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->SlopeEnergy[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->StressTensor[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->InterfacePressure[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Flux_mass[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Flux_diff[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Flux_energy[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForWorkPatch (&(CurrentFluidPatch->Momentum[i]));
  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      ReserveForWorkPatch (&(CurrentFluidPatch->Flux_mom[j][i]));
    }
  }
  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      ReserveForWorkPatch (&(CurrentFluidPatch->Flux_diff_mom[j][i]));
    }
  }

  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      ReserveForWorkPatch (&(CurrentFluidPatch->SlopeVelocity[j][i]));
    }
  }
  if (!Isothermal)
    ReserveForWorkPatch (&(CurrentFluidPatch->Energy_tot));
  ReserveForWorkPatch (&(CurrentFluidPatch->Divergence));
}



void ReserveForSecondaryWorkPatch (ptr)
     real **ptr;
{
  if (!AllocPhase2)
    NField2++;
  else {
    *ptr = FieldStart2+Length2*FieldIndex2;
    FieldIndex2++;
  }
}


void MakeSecondaryFluidPatch ()
{
  long i, j;
  ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Density));
  ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Energy));
  ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Potential));
  ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->SourceRhoPred));
  ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Density_Pred));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Velocity[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Velocity_Pred[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->SourceVelocity[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Accel[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->SlopeDensity[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->SlopeEnergy[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->StressTensor[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->InterfacePressure[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Flux_mass[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Flux_diff[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Flux_energy[i]));
  for (i = 0; i < NDIM; i++)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Momentum[i]));
  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Flux_mom[j][i]));
    }
  }
  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Flux_diff_mom[j][i]));
    }
  }
  for (i = 0; i < NDIM; i++) {
    for (j = 0; j < NDIM; j++) {
      ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->SlopeVelocity[j][i]));
    }
  }
  if (!Isothermal)
    ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Energy_tot));
   ReserveForSecondaryWorkPatch (&(SecondaryFluidPatch->Divergence));
}

void CreateCurrentFluidPatch (fp)
     FluidPatch *fp;
{
  long size[3], Size, i;
   for (i = 0; i < 3; i++)
    size[i] = fp->desc->gncell[i];
  Size = size[0]*size[1]*size[2];
  if (Size == PreviousLength) return;
  if (NeverAllocated) {
    AllocPhase = 0;
    MakeCurrentFluidPatch ();
    AllocPhase = 1;
    CurrentFluidPatch = (FluidWork *)prs_malloc (sizeof(FluidWork));
  }
  if (Size > MaxWorkSize) {
    if (!NeverAllocated)
      free (FieldStart);
    FieldStart = prs_malloc (Size*NField*sizeof(real));
    MaxWorkSize = Size;
    NeverAllocated = FALSE;
  }
  FieldIndex = 0;
  Length = Size;
  PreviousLength = Length;
  MakeCurrentFluidPatch ();
}

void CreateSecondaryFluidPatch (fp)
     FluidPatch *fp;
{
  long size[3], Size, i;
   for (i = 0; i < 3; i++)
    size[i] = fp->desc->gncell[i];
  Size = size[0]*size[1]*size[2];
  if (Size == PreviousLength2) return;
  if (NeverAllocated2) {
    AllocPhase2 = 0;
    MakeSecondaryFluidPatch ();
    AllocPhase2 = 1;
    SecondaryFluidPatch = (FluidWork *)prs_malloc (sizeof(FluidWork));
  }
  if (Size > MaxWorkSize2) {
    if (!NeverAllocated2)
      free (FieldStart2);
    FieldStart2 = prs_malloc (Size*NField2*sizeof(real));
    MaxWorkSize2 = Size;
    NeverAllocated2 = FALSE;
  }
  FieldIndex2 = 0;
  Length2 = Size;
  PreviousLength2 = Length;
  MakeSecondaryFluidPatch ();
}

void ReallocateSpaceForCurrentFluidPatch (fp)
     FluidPatch *fp;
{
  CreateCurrentFluidPatch (fp);
  CurrentFluidPatch->desc = fp->desc;
  CurrentFluidPatch->Fluid = fp;
  CurrentFluidPatch->Fluid->next = fp->next;
  CurrentFluidPatch->next = fp->next;
  CurrentFluidPatch->prev = fp->prev;
}

void ReallocateSpaceForSecondaryFluidPatch (fp)
     FluidPatch *fp;
{
  CreateSecondaryFluidPatch (fp);
  SecondaryFluidPatch->desc = fp->desc;
  SecondaryFluidPatch->Fluid = fp;
  SecondaryFluidPatch->Fluid->next = fp->next;
  SecondaryFluidPatch->next = fp->next;
}

void SendToCurrent (fp)
     FluidPatch *fp;
{
  long gncell[3], stride[3], dm, size;
  ReallocateSpaceForCurrentFluidPatch (fp);
  getgridsize (fp->desc, gncell, stride);
  size = gncell[0]*gncell[1]*gncell[2]*sizeof(real);
  memcpy (CurrentFluidPatch->Density, fp->Density->Field, (size_t)size);
  memcpy (CurrentFluidPatch->Energy, fp->Energy->Field, (size_t)size);
  for (dm = 0; dm < NDIM; dm++) {
    memcpy (CurrentFluidPatch->Velocity[dm], fp->Velocity->Field[dm], (size_t)size);
  }
  if (EXTERNALPOTENTIAL == YES) {
    if ((GlobalDate >= DatePotentialConstant) && (fp->PotentialSet)){
      memcpy(CurrentFluidPatch->Potential, fp->Potential->Field, (size_t)size);
    }
    else {
      ComputeExternalPotential (GlobalDate, fp, 0.0, CurrentFluidPatch->Potential, EVERYWHERE);
    }
  }
  /* Note that the potential in the above expression is evaluated at
     t=0; as such it is a time independent external potential */
}

void SendToSecondary (fp)
     FluidPatch *fp;
{
  long gncell[3], stride[3], dm, size;
  ReallocateSpaceForSecondaryFluidPatch (fp);
  getgridsize (fp->desc, gncell, stride);
  size = gncell[0]*gncell[1]*gncell[2]*sizeof(real);
  memcpy (SecondaryFluidPatch->Density, fp->Density->Field, (size_t)size);
  memcpy (SecondaryFluidPatch->Energy, fp->Energy->Field, (size_t)size);
  for (dm = 0; dm < NDIM; dm++) {
    memcpy (SecondaryFluidPatch->Velocity[dm], fp->Velocity->Field[dm], (size_t)size);
  }
  if (EXTERNALPOTENTIAL == YES) {
    if ((GlobalDate >= DatePotentialConstant) && (fp->PotentialSet)){
      memcpy(SecondaryFluidPatch->Potential, fp->Potential->Field, (size_t)size);
    }
    else {
      ComputeExternalPotential (GlobalDate, fp, 0.0, SecondaryFluidPatch->Potential, EVERYWHERE);
    }
  }
  /* Note that the potential in the above expression is evaluated at
     t=0; as such it is a time independent external potential */
}

void CurrentToPatch (fp)
     FluidPatch *fp;
{
  long gncell[3], stride[3], dm, size;
  ReallocateSpaceForCurrentFluidPatch (fp);
  getgridsize (fp->desc, gncell, stride);
  size = gncell[0]*gncell[1]*gncell[2]*sizeof(real);
  memcpy (fp->Density->Field, CurrentFluidPatch->Density, (size_t)size);
  memcpy (fp->Energy->Field, CurrentFluidPatch->Energy, (size_t)size);
  for (dm = 0; dm < NDIM; dm++) {
    memcpy (fp->Velocity->Field[dm], CurrentFluidPatch->Velocity[dm], (size_t)size);
  }
  if (EXTERNALPOTENTIAL == YES) {
    if (!((GlobalDate >= DatePotentialConstant) && (fp->PotentialSet))) {
      memcpy(fp->Potential->Field, CurrentFluidPatch->Potential, (size_t)size);
      /* This  is  meant to  store  the  potential  field in  the  fluid
	 patch. It will be used e.g. in the torque calculations */
      fp->PotentialSet = YES;
    }
  }
}
