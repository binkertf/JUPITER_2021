#include "jupiter.h"

extern long TimeStepRatio[REFINEUPPERLIMIT];
void CorrectFluxesFromFinerLevel ()
{
  long nvar, gncell[3], stride[3], i, dim, side, dim1, dim2;
  long imin[3], imax[3], size[3], j, k, m, l, n;
  long pres_varindex = 0;
  FluidWork *fw;
  Communicator *com;
  real *dest[8][3], *buffer, *inter[3];
  nvar = 2+NDIM;		/* Density (or mass flux) + velocities (or momentum flux) */
  if (!Isothermal) nvar+=2;	/* Energy of energy flux */
  fw = CurrentFluidPatch;
  com = ComListFlux;
  getgridsize (fw->desc, gncell, stride);
  for (i = 0; i < 3; i++)
    inter[i] = fw->desc->InterSurface[i];
  for (j = 0; j < NDIM; j++) {
    dest[0][j] = fw->Flux_mass[j];
    for (i = 1; i <= NDIM; i++)
      dest[i][j] = fw->Flux_mom[i-1][j];
    pres_varindex = NDIM+1;
    dest[pres_varindex][j] = fw->InterfacePressure[j];
    if (!Isothermal) {
      dest[NDIM+2][j] = fw->Flux_energy[j];
      dest[NDIM+3][j] = fw->Flux_tot_energy[j];
    }
  }
  while (com != NULL) {
    if (com->destg == fw->desc) {
      dim = com->facedim;
      side = com->faceside;
      dim1 = (dim == 0);
      dim2 = 2 - (dim == 2);
      for (i = 0; i < 3; i++) {
	imin[i] = com->imin_dest[i];
	imax[i] = com->imax_dest[i];
	size[i] = imax[i]-imin[i];
      }
      buffer = (com->buffer)+(fw->Fluid->FluidRank)*nvar*(com->size);
      i = (side == INF ? imax[dim] : imin[dim]);
      for (j = imin[dim1]; j < imax[dim1]; j++) {
	for (k = imin[dim2]; k < imax[dim2]; k++) {
	  m=i*stride[dim]+j*stride[dim1]+k*stride[dim2];
	  for (l = 0; l < nvar; l++) {
	    n = (j-imin[dim1])+(k-imin[dim2])*size[dim1];
	      dest[l][dim][m] = buffer[n+l*size[dim1]*size[dim2]];
	    if (l == pres_varindex)
	      dest[l][dim][m] /= (real)(TimeStepRatio[fw->desc->level])*inter[dim][m];
	  }
	}
      }
    }
    com = com->next;
  }
}

