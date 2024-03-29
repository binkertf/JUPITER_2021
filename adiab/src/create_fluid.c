#include "jupiter.h"

ScalarField *CreateScalarField (desc, name, ptr)
     tGrid_CPU *desc;
     char *name;
     real *ptr;
{
  char *Name;
  char string[MAXLINELENGTH];
  ScalarField *sf;
  long size[3], i;
  real *Field;
  sf = prs_malloc (sizeof(ScalarField));
  if (strlen(name) > MAXNAMELENGTH) {
    sprintf (string, "%s is longer than the %ld characters limit.", name, MAXNAMELENGTH);
    prs_error (string);
  }
  Name = prs_malloc(sizeof(char)*(MAXNAMELENGTH+1));
  strcpy (Name, name);
  sf->Name = Name;
  sf->desc = desc;
  for (i = 0; i < 3; i++)
    size[i] = desc->gncell[i];
  Field = ptr;
  sf->Field = Field;
  return sf;
}

void FreeScalarField (sf)
     ScalarField *sf;
{
  free (sf->Name);
  free (sf);
}

VectorField *CreateVectorField (desc, name, ptr)
     tGrid_CPU *desc;
     char *name;
     real *ptr;
{
  char *Name;
  char string[MAXLINELENGTH];
  VectorField *vf;
  long size[3], i;
  real *Field;
  vf = prs_malloc (sizeof(VectorField));
  if (strlen(name) > MAXNAMELENGTH) {
    sprintf (string, "%s is longer than the %ld characters limit.", name, MAXNAMELENGTH);
    prs_error (string);
  }
  Name = prs_malloc(sizeof(char)*(MAXNAMELENGTH+1));
  strcpy (Name, name);
  vf->Name = Name;
  vf->desc = desc;
  for (i = 0; i < 3; i++)
    size[i] = desc->gncell[i];
  Field = ptr;
  for (i = 0; i < NDIM; i++) {
    vf->Field[i] = Field+i*size[0]*size[1]*size[2]; /* The velocity fields are contiguous in memory */
   }
  return vf;
}

void FreeVectorField (vf)
     VectorField *vf;
{
  free (vf->Name);
  free (vf);
}

InterfaceFlux *CreateInterfaceFlux (desc)
     tGrid_CPU *desc;
{
  InterfaceFlux *iflux;
  real *field[3][2];
  long dim, side, dim1, dim2, size[3];
  iflux = (InterfaceFlux *)prs_malloc (sizeof(InterfaceFlux));
  iflux->desc = desc;
  for (dim = 0; dim < 3; dim++)	/* Do not change this 3 (value 1 is
				   needed below for extra
				   dimensions) */
    size[dim] = desc->gncell[dim];
  for (dim = 0; dim < NDIM; dim++) {
    dim1 = (dim == 0);
    dim2 = 2 - (dim == 2);
    for (side = INF; side <= SUP; side++) {
      field[dim][side] = prs_malloc(sizeof(real)*size[dim1]*size[dim2]);
      iflux->Flux[dim][side] = field[dim][side];
    }
  }
  return iflux;
}

void FreeInterfaceFlux (iflux)
     InterfaceFlux *iflux;
{
  long dim, side;
  for (dim = 0; dim < NDIM; dim++) {
    for (side = INF; side <= SUP; side++) {
      free (iflux->Flux[dim][side]);
    }
  }
  free (iflux);
}

PressureFaces *CreatePressureFaces (desc)
     tGrid_CPU *desc;
{
  PressureFaces *pres;
  real *field[3][2];
  long dim, side, dim1, dim2, size[3];
  pres = (PressureFaces *)prs_malloc (sizeof(PressureFaces));
  pres->desc = desc;
  for (dim = 0; dim < 3; dim++)	/* Do not change this 3 (value 1 is
				   needed below for extra
				   dimensions) */
    size[dim] = desc->gncell[dim];
  for (dim = 0; dim < NDIM; dim++) {
    dim1 = (dim == 0);
    dim2 = 2 - (dim == 2);
    for (side = INF; side <= SUP; side++) {
      field[dim][side] = prs_malloc(sizeof(real)*size[dim1]*size[dim2]);
      pres->Pressure[dim][side] = field[dim][side];
    }
  }
  return pres;
}

void FreePressureFaces (pres)
     PressureFaces *pres;
{
  long dim, side;
  for (dim = 0; dim < NDIM; dim++) {
    for (side = INF; side <= SUP; side++) {
      free (pres->Pressure[dim][side]);
    }
  }
  free (pres);
}

FluidPatch *CreateFluidPatch (desc, name, initcode, initcodeeq)
     tGrid_CPU *desc;
     char *name, *initcode, *initcodeeq;
{
  char *Name;
  real *StartField;
  long dim, i, size[3], Size;
  long MeridianSize;
  ScalarField *Density, *Energy, *TotalEnergy, *Potential, *OpticalDepth, *EnergyRad;
  ScalarField *StellarHeating, *EradDeriv, *Temperature, *Opacity, *TauCell, *Gamma;
  VectorField *Velocity;
  InterfaceFlux *MassFlux, *MomentumFlux[3], *EnergyFlux, *TotalEnergyFlux, *DiffFlux;
  PressureFaces *Pressure;
  FluidPatch *patch;
  int n_stellar; // number of additional properties that are called if Stellar == True

  patch = prs_malloc (sizeof(FluidPatch));
  patch->desc = desc;
  if (strlen(name) > MAXNAMELENGTH)
    prs_error ("%s is longer than the %ld characters limit.", name, MAXNAMELENGTH);
  Name = prs_malloc(sizeof(char)*(MAXNAMELENGTH+1));
  strcpy (Name, name);
  patch->Name = Name;
  for (i = 0; i < 3; i++)
    size[i] = desc->gncell[i];
  Size = size[0] * size[1] * size[2];
  n_stellar = 8;   
  StartField = prs_malloc (sizeof(real)*Size*(1+1+NDIM+1+1+n_stellar*(Stellar == YES ? 1 : 0))); // number of fields below
  /* Global contiguous allocation */
  patch->StartField = StartField;
  Density = CreateScalarField (desc, "density", StartField);
  /* Density MUST be the first field (see function ResetPatch below) */
  patch->Ptr[_Density_] = Density->Field;
  Energy = CreateScalarField (desc, "energy", StartField+Size);
  patch->Ptr[_Energy_] = Energy->Field;
  Velocity = CreateVectorField (desc, "velocity", StartField+2*Size);
  Potential = CreateScalarField (desc, "potential", StartField+(2+NDIM)*Size);
  TotalEnergy = CreateScalarField (desc, "totalenergy", StartField+(3+NDIM)*Size);
  if (Stellar) {
    OpticalDepth   = CreateScalarField (desc, "tau", StartField+(4+NDIM)*Size);
    EnergyRad      = CreateScalarField (desc, "erad", StartField+(5+NDIM)*Size);
    StellarHeating = CreateScalarField (desc, "stheat", StartField+(6+NDIM)*Size);
    EradDeriv      = CreateScalarField (desc, "erderiv", StartField+(7+NDIM)*Size);
    Temperature    = CreateScalarField (desc, "temperature", StartField+(8+NDIM)*Size);
    Opacity        = CreateScalarField (desc, "opacity", StartField+(9+NDIM)*Size);
    TauCell        = CreateScalarField (desc, "taucell", StartField+(10+NDIM)*Size);
    Gamma          = CreateScalarField (desc, "gamma", StartField+(11+NDIM)*Size);
  }
  patch->Ptr[_Tot_Energy_] = TotalEnergy->Field;
  patch->Ptr[_Potential_] = Potential->Field;
  patch->Ptr[_Velocity_] = Velocity->Field[0];
  patch->Ptr[_Velocity_+1] = Velocity->Field[1];
  patch->Ptr[_Velocity_+2] = Velocity->Field[2];
  if (Stellar) {
    patch->Ptr[_Tau_] =  OpticalDepth->Field;
    patch->Ptr[_Erad_] = EnergyRad->Field;
    patch->Ptr[_EradDeriv_] = EradDeriv->Field;
    patch->Ptr[_Temperature_] = Temperature->Field;
    patch->Ptr[_Opacity_] = Opacity->Field;
    patch->Ptr[_TauCell_] = TauCell->Field;
    patch->Ptr[_Gamma_] = Gamma->Field; 
  }
  MassFlux = CreateInterfaceFlux (desc);
  EnergyFlux = CreateInterfaceFlux (desc);
  TotalEnergyFlux = CreateInterfaceFlux (desc);
  DiffFlux = CreateInterfaceFlux (desc);
  if (KEPLERIAN) {		/* We create 3 contiguous arrays that
				   control the azimuthal flux
				   corrections for submeshes */
    MeridianSize = size[_RAD_]*size[_COLAT_];
    patch->MassFluxCorrection1     = prs_malloc (sizeof(real)*MeridianSize*6L);
    patch->MassFluxCorrection2     = patch->MassFluxCorrection1+1L*MeridianSize;
    patch->MomentumFluxCorrection1 = patch->MassFluxCorrection1+2L*MeridianSize;
    patch->MomentumFluxCorrection2 = patch->MassFluxCorrection1+3L*MeridianSize;
    patch->EnergyFluxCorrection1   = patch->MassFluxCorrection1+4L*MeridianSize;
    patch->EnergyFluxCorrection2   = patch->MassFluxCorrection1+5L*MeridianSize;
  }
  for (dim = 0; dim < NDIM; dim++)
    MomentumFlux[dim] = CreateInterfaceFlux (desc);
  Pressure = CreatePressureFaces (desc);
  patch->Pressure = Pressure;
  patch->Density = Density;
  patch->Energy = Energy;
  patch->TotalEnergy = TotalEnergy;
  patch->Velocity = Velocity;
  if (Stellar) {
    patch->OpticalDepth = OpticalDepth;
    patch->EnergyRad    = EnergyRad;
    patch->StellarHeating = StellarHeating;
    patch->EradDeriv = EradDeriv;
    patch->Temperature = Temperature;
    patch->Opacity = Opacity;
    patch->TauCell = TauCell;
    patch->Gamma = Gamma;
  }
  patch->MassFlux = MassFlux;
  patch->DiffFlux = DiffFlux;
  patch->EnergyFlux = EnergyFlux;
  patch->TotalEnergyFlux = TotalEnergyFlux;
  patch->PreviousEradExists = NO;
  for (dim = 0; dim < NDIM; dim++)
    patch->MomentumFlux[dim] = MomentumFlux[dim];
  patch->Potential = Potential;
  patch->PotentialSet = FALSE;
  findcodeinit (initcode, &(patch->InitCode));
  findcodeinit (initcodeeq, &(patch->InitCodeEq));
  return patch;
}

void FreeFluidPatch (patch)
     FluidPatch *patch;
{
  long dim;
  free (patch->Name);
  FreeScalarField (patch->Density);
  FreeScalarField (patch->Energy);
  FreeScalarField (patch->TotalEnergy);
  FreeVectorField (patch->Velocity);
  if (Stellar) {
    FreeScalarField (patch->OpticalDepth);
    FreeScalarField (patch->EnergyRad);
    FreeScalarField (patch->StellarHeating);
    FreeScalarField (patch->Temperature);
    FreeScalarField (patch->Opacity);
    FreeScalarField (patch->TauCell);
    FreeScalarField (patch->Gamma);
  }
  FreeInterfaceFlux (patch->MassFlux);
  FreeInterfaceFlux (patch->DiffFlux);
  FreeInterfaceFlux (patch->EnergyFlux);
  FreeInterfaceFlux (patch->TotalEnergyFlux);
  FreePressureFaces (patch->Pressure);
  for (dim = 0; dim < NDIM; dim++)
    FreeInterfaceFlux (patch->MomentumFlux[dim]);
  FreeScalarField (patch->Potential);
  free (patch->StartField);
  free (patch);
}

void ResetPatch (patch)
     FluidPatch *patch;
{
  long i, Size=1;
  long nb;
  nb=2+NDIM;
  if (!Isothermal) nb = 3+NDIM;
  if (Stellar) nb+=7;
  for (i = 0; i < 3; i++)
    Size *= patch->desc->gncell[i];
  SetToZero (patch->Density->Field, nb*Size);
  /* The above '2+NDIM'/'3+NDIM' is because all fields are contiguous in memory
     and 'Density' is the first one */
}
