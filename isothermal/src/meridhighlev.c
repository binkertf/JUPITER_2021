/** \file meridhighlev.c

 Find the highest level of the mesh system on a given ring. Find the
highest level of the mesh system within the ring that intersects the
coarsest mesh on the underlying coarsest zone. This is meant to
perform integrals for the azimuthal fluxes corrections with an
appropriate number of integration substeps (see meridcorrflux.c). For
a given fluid patch, the corresponding data set is a 2D map that has
the size of a meridian cut of the mesh ("Nrad" * "Ncolat"), and which is
stored in the field "MaxLevMeridianProj". Note that communications between
processing elements are required in order to get the highest level over
all processing elements (hence the use of the arrays "CoarseMeridian"
and "CoarseMeridianG" -G like global-). 

 */

#include "jupiter.h"

void GetMaxLevProjectedMeridian () {  
  tGrid_CPU *g;
  long gncell[3], stride[3], n1, n2;
  long m, i, j, *ml, levg, ic, jc, mc;
  long k, Incr0[3], Incrl[3], ia, ja;
  int *CoarseMeridian, *CoarseMeridianG;
  long nradc, ncolc;
  nradc = Ncell0[_RAD_]+2*Nghost[_RAD_];
  ncolc = Ncell0[_COLAT_]+2*Nghost[_COLAT_];
  CoarseMeridian = prs_malloc (sizeof(int)*nradc*ncolc);
  CoarseMeridianG = prs_malloc (sizeof(int)*nradc*ncolc);
  for (i = 0; i < nradc * ncolc; i++) {
    CoarseMeridian[i] = 0;
    CoarseMeridianG[i] = 0;
  }
  /* Note that the following only works for grids that do not exceed
     the template */
  for (k = 0; k < 3; k++)
    Incr0[k] = (ncorner_max0[k]-ncorner_min0[k])/Ncell0[k];
  g = Grid_CPU_list;
  while (g != NULL) {		/* For each grid in the current process */
    if (g->cpu == CPU_Rank) {
     getgridsize (g, gncell, stride);
     levg = g->level;
      for (k = 0; k < 3; k++)
	Incrl[k] = (g->gncorner_max[k]-g->gncorner_min[k])/gncell[k];
      for (i = 0; i < gncell[_RAD_]; i++)  {
	for (j = 0; j < gncell[_COLAT_]; j++)  {
	  /* Where are we on the absolute grid ? */
	  ia = g->gncorner_min[_RAD_]+i*Incrl[_RAD_];
	  ja = g->gncorner_min[_COLAT_]+j*Incrl[_COLAT_];
	  /* And now on the coarsest grid ?*/
	  ia -= ncorner_min0[_RAD_];
	  ja -= ncorner_min0[_COLAT_];
	  ic = ia/Incr0[_RAD_];	/* Improve this for negative values */
	  jc = ja/Incr0[_COLAT_]; /* idem  */
	  ic += Nghost[_RAD_];
	  jc += Nghost[_COLAT_];
	  mc = ic*ncolc+jc;
	  if (CoarseMeridian[mc] < (int)levg)
	    CoarseMeridian[mc] = (int)levg;
	}
      }
    }
    g = g->next;
  }
  MPI_Allreduce (CoarseMeridian, CoarseMeridianG, (int)(nradc*ncolc), MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  g = Grid_CPU_list;
  while (g != NULL) {		/* For each grid in the current process */
    if (g->cpu == CPU_Rank) {
      getgridsize (g, gncell, stride);
      ml = g->MaxLevMeridianProj;
      n1 = gncell[_RAD_];
      n2 = gncell[_COLAT_];
      for (k = 0; k < 3; k++)
	Incrl[k] = (g->gncorner_max[k]-g->gncorner_min[k])/gncell[k];
      for (i = 0; i < n1; i++)  {
	for (j = 0; j < n2; j++)  {
	  /* Where are we on the absolute grid ? */
	  ia = g->gncorner_min[_RAD_]+i*Incrl[_RAD_];
	  ja = g->gncorner_min[_COLAT_]+j*Incrl[_COLAT_];
	  /* And now on the coarsest grid ?*/
	  ia -= ncorner_min0[_RAD_];
	  ja -= ncorner_min0[_COLAT_];
	  ic = ia/Incr0[_RAD_];
	  jc = ja/Incr0[_COLAT_];
	  ic += Nghost[_RAD_];
	  jc += Nghost[_COLAT_];
	  mc = ic*ncolc+jc;
	  m  = i+j*n1;
	  if (_COLAT_ < _RAD_) m =j+i*n2;
	  ml[m] = (long)CoarseMeridianG[mc];
	}
      }
    }
    g = g->next;
  }
  free (CoarseMeridian);
  free (CoarseMeridianG);
}
