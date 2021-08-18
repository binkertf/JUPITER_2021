/** \file hstatwrite.c

Output the squeezed arrays that contain information about the
equilibrium state. The squeezed arrays are unfold so as to have
dimensions matching those of the normal arrays, so that they can be
read straightforwardly by the standard IDL routines.

 */

#include "jupiter.h"

void WriteSqueezedField (Fluid, Field, number)
     FluidPatch *Fluid;
     SqueezedField *Field;
     long number;
{
  FILE *hdl;
  long gncell[3], Size, i, j, k, ms, m;
  long strides[3], stride[3], ex[3], ngh[3];
  real *buffer;
  char filename[MAXLINELENGTH];
  sprintf (filename, "%s%s%ld_%ld_%ld.dat", Fluid->Name, Field->Name,	\
	   number, Field->desc->number, Field->desc->level);
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = (WriteGhosts ? Nghost[i] : 0);
    ex[i] = Nghost[i]-ngh[i];
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
    strides[i]= Field->stride[i];
  }
  Size = (gncell[0]-2*ex[0])*(gncell[1]-2*ex[1])*(gncell[2]-2*ex[2]);
  m = 0;
  buffer = prs_malloc (Size*sizeof(real));
  for (k = ex[2]; k < gncell[2]-ex[2]; k++) {
    for (j = ex[1]; j < gncell[1]-ex[1]; j++) {
      for (i = ex[0]; i < gncell[0]-ex[0]; i++) {
	ms = i*strides[0]+j*strides[1]+k*strides[2];
	buffer[m] = Field->field[ms];
	m++;
      }
    }
  }
  fwrite (buffer, sizeof(real), Size, hdl);
  free (buffer);
  fclose (hdl);
}

void WriteHydroStat (number)
     long number;
{
  tGrid_CPU *item;
  FluidPatch *fluid;
  long i;
  item = Grid_CPU_list;
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      fluid = item->Fluid;
      while (fluid !=NULL) {
	WriteSqueezedField (fluid, fluid->Rho_eq_c, number);
	for (i = 0; i < 3; i++) {
	  if (CorrHydroStat[i]) {
	    WriteSqueezedField (fluid, fluid->Rho_eq_i[i], number);
	    WriteSqueezedField (fluid, fluid->Pres_eq_i[i], number);
	    WriteSqueezedField (fluid, fluid->Source_corr[i], number);
	    WriteSqueezedField (fluid, fluid->Cs2_i[i], number);
	  }
	}
	fluid = fluid->next;
      }
    }
    item = item->next;
  }
}
