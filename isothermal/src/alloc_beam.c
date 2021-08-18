#include "jupiter.h"

static long MaxBeamLength=0;
static boolean NeverAllocated = TRUE;
static long AllocPhase;
static long Length;
static long NVar, VarIndex;
static real *BeamStart;
static long PreviousLength=0;

void ReserveForBeam_PLM (ptr)
     real **ptr;
{
  if (!AllocPhase)
    NVar++;
  else {
    *ptr = BeamStart+Length*VarIndex;
    VarIndex++;
  }
}

void MakeBeam_PLM (beam)
Beam *beam;
{
  long i;
  ReserveForBeam_PLM (&(beam->rho));
  ReserveForBeam_PLM (&(beam->rho_pred));
  ReserveForBeam_PLM (&(beam->u));
  ReserveForBeam_PLM (&(beam->u_pred));
  ReserveForBeam_PLM (&(beam->cs));
  for (i = 0; i < 3; i++) {	/* Could be NDIM? */
    ReserveForBeam_PLM (&(beam->slope_rho[i]));
    ReserveForBeam_PLM (&(beam->slope_u[i]));
    ReserveForBeam_PLM (&(beam->slope_u[i]));
    ReserveForBeam_PLM (&(beam->slope_v_perp[0][i]));
    ReserveForBeam_PLM (&(beam->slope_v_perp[1][i]));
    ReserveForBeam_PLM (&(beam->momentum_flux[i]));
    ReserveForBeam_PLM (&(beam->momentum_diff_flux[i]));
    if (!Isothermal)
      ReserveForBeam_PLM (&(beam->slope_energy[i]));
  }
  if (!Isothermal)
    ReserveForBeam_PLM (&beam->energy_flux);
  ReserveForBeam_PLM (&(beam->center));
  ReserveForBeam_PLM (&(beam->rawcoord));
  ReserveForBeam_PLM (&(beam->edge));
  ReserveForBeam_PLM (&(beam->srcrho));
  ReserveForBeam_PLM (&(beam->intersurface));
  ReserveForBeam_PLM (&(beam->InvVolume));
  ReserveForBeam_PLM (&(beam->mass_flux));
  ReserveForBeam_PLM (&(beam->diff_flux));
  ReserveForBeam_PLM (&(beam->rhoL));
  ReserveForBeam_PLM (&(beam->rhoR));
  ReserveForBeam_PLM (&(beam->uL));
  ReserveForBeam_PLM (&(beam->uR));
  ReserveForBeam_PLM (&(beam->eL));
  ReserveForBeam_PLM (&(beam->eR));
  for (i = 0; i < 2; i++) {
    ReserveForBeam_PLM (&(beam->v_perp[i]));
    ReserveForBeam_PLM (&(beam->v_perp_pred[i]));
    ReserveForBeam_PLM (&(beam->v_perp_L[i]));
    ReserveForBeam_PLM (&(beam->v_perp_R[i]));
    ReserveForBeam_PLM (&(beam->metperp[i]));
    ReserveForBeam_PLM (&(beam->source_perp[i]));
    ReserveForBeam_PLM (&(beam->MomCorrRatioL[i]));
    ReserveForBeam_PLM (&(beam->MomCorrRatioR[i]));
  }
  ReserveForBeam_PLM (&(beam->pressure_godunov));
  ReserveForBeam_PLM (&(beam->source));
  ReserveForBeam_PLM (&(beam->HS_cent_rho));
  ReserveForBeam_PLM (&(beam->HS_int_rho));
  ReserveForBeam_PLM (&(beam->cs2i));
}

void AllocBeam_PLM (beam, length)
     Beam *beam;
     long length;
{
  if (length == PreviousLength) return;
  if (NeverAllocated) {		/* We need to know how many */
    /* variables  the  beam will  contain.  This  number of  variables
     cannot therefore depend upon the coordinate */
    AllocPhase = 0;
    MakeBeam_PLM (beam);
    AllocPhase = 1;
  }
  if (length > MaxBeamLength) {
    if (!NeverAllocated)
      free (BeamStart);
    BeamStart = prs_malloc (length*NVar*sizeof(real));
    MaxBeamLength = length;
    NeverAllocated = FALSE;
    pInfo ("Beam contains %ld variables and each field has length at most %ld\n", NVar, length);
  }
  VarIndex = 0;
  Length = length;
  PreviousLength = Length;
  MakeBeam_PLM (beam);
}


static long MaxBeamLength2=0;
static boolean NeverAllocated2 = TRUE;
static long AllocPhase2;
static long Length2;
static long NVar2, VarIndex2;
static real *BeamStart2;
static long PreviousLength2=0;

void ReserveForSecondaryBeam_PLM (ptr)
     real **ptr;
{
  if (!AllocPhase2)
    NVar2++;
  else {
    *ptr = BeamStart2+Length2*VarIndex2;
    VarIndex2++;
  }
}

void MakeSecondaryBeam_PLM (beam)
Beam *beam;
{
  long i;
  ReserveForSecondaryBeam_PLM (&(beam->rho));
  ReserveForSecondaryBeam_PLM (&(beam->rho_pred));
  ReserveForSecondaryBeam_PLM (&(beam->u));
  ReserveForSecondaryBeam_PLM (&(beam->u_pred));
  ReserveForSecondaryBeam_PLM (&(beam->cs));
  for (i = 0; i < 3; i++) {	/* Could be NDIM? */
    ReserveForSecondaryBeam_PLM (&(beam->slope_rho[i]));
    ReserveForSecondaryBeam_PLM (&(beam->slope_u[i]));
    ReserveForSecondaryBeam_PLM (&(beam->slope_u[i]));
    ReserveForSecondaryBeam_PLM (&(beam->slope_v_perp[0][i]));
    ReserveForSecondaryBeam_PLM (&(beam->slope_v_perp[1][i]));
    ReserveForSecondaryBeam_PLM (&(beam->momentum_flux[i]));
    ReserveForSecondaryBeam_PLM (&(beam->momentum_diff_flux[i]));
    if (!Isothermal)
      ReserveForSecondaryBeam_PLM (&(beam->slope_energy[i]));
  }
  if (!Isothermal)
    ReserveForSecondaryBeam_PLM (&beam->energy_flux);
  ReserveForSecondaryBeam_PLM (&(beam->center));
  ReserveForSecondaryBeam_PLM (&(beam->rawcoord));
  ReserveForSecondaryBeam_PLM (&(beam->edge));
  ReserveForSecondaryBeam_PLM (&(beam->srcrho));
  ReserveForSecondaryBeam_PLM (&(beam->intersurface));
  ReserveForSecondaryBeam_PLM (&(beam->InvVolume));
  ReserveForSecondaryBeam_PLM (&(beam->mass_flux));
  ReserveForSecondaryBeam_PLM (&(beam->diff_flux));
  ReserveForSecondaryBeam_PLM (&(beam->rhoL));
  ReserveForSecondaryBeam_PLM (&(beam->rhoR));
  ReserveForSecondaryBeam_PLM (&(beam->uL));
  ReserveForSecondaryBeam_PLM (&(beam->uR));
  ReserveForSecondaryBeam_PLM (&(beam->eL));
  ReserveForSecondaryBeam_PLM (&(beam->eR));
  for (i = 0; i < 2; i++) {
    ReserveForSecondaryBeam_PLM (&(beam->v_perp[i]));
    ReserveForSecondaryBeam_PLM (&(beam->v_perp_pred[i]));
    ReserveForSecondaryBeam_PLM (&(beam->v_perp_L[i]));
    ReserveForSecondaryBeam_PLM (&(beam->v_perp_R[i]));
    ReserveForSecondaryBeam_PLM (&(beam->metperp[i]));
    ReserveForSecondaryBeam_PLM (&(beam->source_perp[i]));
    ReserveForSecondaryBeam_PLM (&(beam->MomCorrRatioL[i]));
    ReserveForSecondaryBeam_PLM (&(beam->MomCorrRatioR[i]));
  }
  ReserveForSecondaryBeam_PLM (&(beam->pressure_godunov));
  ReserveForSecondaryBeam_PLM (&(beam->source));
  ReserveForSecondaryBeam_PLM (&(beam->HS_cent_rho));
  ReserveForSecondaryBeam_PLM (&(beam->HS_int_rho));
  ReserveForSecondaryBeam_PLM (&(beam->cs2i));
}

void AllocSecondaryBeam_PLM (beam, length)
     Beam *beam;
     long length;
{
  if (length == PreviousLength2) return;
  if (NeverAllocated2) {		/* We need to know how many */
    /* variables  the  beam will  contain.  This  number of  variables
     cannot therefore depend upon the coordinate */
    AllocPhase2 = 0;
    MakeSecondaryBeam_PLM (beam);
    AllocPhase2 = 1;
  }
  if (length > MaxBeamLength2) {
    if (!NeverAllocated2)
      free (BeamStart2);
    BeamStart2 = prs_malloc (length*NVar*sizeof(real));
    MaxBeamLength2 = length;
    NeverAllocated2 = FALSE;
    pInfo ("Beam contains %ld variables and each field has length at most %ld\n", NVar2, length);
  }
  VarIndex2 = 0;
  Length2 = length;
  PreviousLength2 = Length2;
  MakeSecondaryBeam_PLM (beam);
}

//-----------------------------------------------------------------------------------------------------------------------------

static long MaxBeamLength3=0;
static boolean NeverAllocated3 = TRUE;
static long AllocPhase3;
static long Length3;
static long NVar3, VarIndex3;
static real *BeamStart3;
static long PreviousLength3=0;

void ReserveForTertiaryBeam_PLM (ptr)
     real **ptr;
{
  if (!AllocPhase3)
    NVar3++;
  else {
    *ptr = BeamStart3+Length3*VarIndex3;
    VarIndex3++;
  }
}

void MakeTertiaryBeam_PLM (beam)
Beam *beam;
{
  long i;
  ReserveForTertiaryBeam_PLM (&(beam->rho));
  ReserveForTertiaryBeam_PLM (&(beam->rho_pred));
  ReserveForTertiaryBeam_PLM (&(beam->u));
  ReserveForTertiaryBeam_PLM (&(beam->u_pred));
  ReserveForTertiaryBeam_PLM (&(beam->cs));
  for (i = 0; i < 3; i++) {	/* Could be NDIM? */
    ReserveForTertiaryBeam_PLM (&(beam->slope_rho[i]));
    ReserveForTertiaryBeam_PLM (&(beam->slope_u[i]));
    ReserveForTertiaryBeam_PLM (&(beam->slope_u[i]));
    ReserveForTertiaryBeam_PLM (&(beam->slope_v_perp[0][i]));
    ReserveForTertiaryBeam_PLM (&(beam->slope_v_perp[1][i]));
    ReserveForTertiaryBeam_PLM (&(beam->momentum_flux[i]));
    ReserveForTertiaryBeam_PLM (&(beam->momentum_diff_flux[i]));
    if (!Isothermal)
      ReserveForTertiaryBeam_PLM (&(beam->slope_energy[i]));
  }
  if (!Isothermal)
    ReserveForTertiaryBeam_PLM (&beam->energy_flux);
  ReserveForTertiaryBeam_PLM (&(beam->center));
  ReserveForTertiaryBeam_PLM (&(beam->rawcoord));
  ReserveForTertiaryBeam_PLM (&(beam->edge));
  ReserveForTertiaryBeam_PLM (&(beam->srcrho));
  ReserveForTertiaryBeam_PLM (&(beam->intersurface));
  ReserveForTertiaryBeam_PLM (&(beam->InvVolume));
  ReserveForTertiaryBeam_PLM (&(beam->mass_flux));
  ReserveForTertiaryBeam_PLM (&(beam->diff_flux));
  ReserveForTertiaryBeam_PLM (&(beam->rhoL));
  ReserveForTertiaryBeam_PLM (&(beam->rhoR));
  ReserveForTertiaryBeam_PLM (&(beam->uL));
  ReserveForTertiaryBeam_PLM (&(beam->uR));
  ReserveForTertiaryBeam_PLM (&(beam->eL));
  ReserveForTertiaryBeam_PLM (&(beam->eR));
  for (i = 0; i < 2; i++) {
    ReserveForTertiaryBeam_PLM (&(beam->v_perp[i]));
    ReserveForTertiaryBeam_PLM (&(beam->v_perp_pred[i]));
    ReserveForTertiaryBeam_PLM (&(beam->v_perp_L[i]));
    ReserveForTertiaryBeam_PLM (&(beam->v_perp_R[i]));
    ReserveForTertiaryBeam_PLM (&(beam->metperp[i]));
    ReserveForTertiaryBeam_PLM (&(beam->source_perp[i]));
    ReserveForTertiaryBeam_PLM (&(beam->MomCorrRatioL[i]));
    ReserveForTertiaryBeam_PLM (&(beam->MomCorrRatioR[i]));
  }
  ReserveForTertiaryBeam_PLM (&(beam->pressure_godunov));
  ReserveForTertiaryBeam_PLM (&(beam->source));
  ReserveForTertiaryBeam_PLM (&(beam->HS_cent_rho));
  ReserveForTertiaryBeam_PLM (&(beam->HS_int_rho));
  ReserveForTertiaryBeam_PLM (&(beam->cs2i));
}

void AllocTertiaryBeam_PLM (beam, length)
     Beam *beam;
     long length;
{
  if (length == PreviousLength3) return;
  if (NeverAllocated3) {		/* We need to know how many */
    /* variables  the  beam will  contain.  This  number of  variables
     cannot therefore depend upon the coordinate */
    AllocPhase3 = 0;
    MakeTertiaryBeam_PLM (beam);
    AllocPhase3 = 1;
  }
  if (length > MaxBeamLength3) {
    if (!NeverAllocated3)
      free (BeamStart3);
    BeamStart3 = prs_malloc (length*NVar*sizeof(real));
    MaxBeamLength3 = length;
    NeverAllocated3 = FALSE;
    pInfo ("Beam contains %ld variables and each field has length at most %ld\n", NVar3, length);
  }
  VarIndex3 = 0;
  Length3 = length;
  PreviousLength3 = Length3;
  MakeTertiaryBeam_PLM (beam);
}
