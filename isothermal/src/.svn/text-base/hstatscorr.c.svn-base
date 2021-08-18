/** \file hstatscorr.c

  Tabulate the (additive) correction to the source term. This
correction is needed to get a strict numerical equilibrium in the
equilibrium configuration.  In addition, some statistics (which may
help to diagnose an incorrect equilibrium) are written in info.log

*/

#include "jupiter.h"

extern char *CoordNames[];

void CorrectSourceTerm (fluid)
     FluidPatch *fluid;
{
  SqueezedField *sc[3], *pi[3];
  real *p_int;
  real *sv[3], *inter[3], *invvol, *rho;
  real corrmean, corrmin, corrmax, corr;
  real rcorrmean, rcorrmin, rcorrmax; /* r for reduced mesh */
  long gncell[3], stride[3], strides[3], foo[3], nbcorr=0;
  long i, j, k, l, idx[3], ms, m, lp1, lp2, rnbcorr=0;
  FluidWork *fw;
  SendToCurrent (fluid);
  fw = CurrentFluidPatch;
  FillSources (PREDICT, EVERYWHERE);
  invvol = fluid->desc->InvVolume;
  rho = fluid->Density->Field;
  getgridsizes (fluid->desc,  foo, stride, gncell, strides);
  for (i = 0; i < 3; i++) {
    sv[i] = fw->SourceVelocity[i];
    pi[i] = fluid->Pres_eq_i[i];
    sc[i] = fluid->Source_corr[i];
    inter[i] = fluid->desc->InterSurface[i];
    if (gncell[i] > 1) gncell[i]--;
  }
  /* In what follows, we make a loop on the dimension and seek the
     dimensions along which hydrostatic equilibrium is enforced. The
     choice of the loop limits is tricky: we want to scan the squeezed
     arrays, not the full arrays, for speed reason. So potentially
     gncell[] can be 1. However, the source velocity arrays, for the
     full grid, is defined on the main grid + 1 ghost only, not two
     ghosts. Therefore, the loop extends over the full arrays, ghosts
     included, and the index which correspond to invariant dimensions
     are shifted by NGH. Later on, the statistics on the source term
     corrections are performed only for the zones that are in the
     active part (+1 GHOST) (provided the corresponding dimension is
     not invariant, ie their gncell is not 1) */
  for (l = 0; l < 3; l++) {
    lp1 = (l == 0);
    lp2 = 2 - (l == 2);
    if (CorrHydroStat[l]) {
      rcorrmean = corrmean = rcorrmax = corrmax = 0.0;
      rcorrmin = corrmin = 1e30;
      rnbcorr = nbcorr = 0;
      p_int = pi[l]->field;
      for (k = 0; k < gncell[2]; k++) {
	idx[2] = k;
	if (SqueezeDim[2]) idx[2] += NGH;
	for (j = 0; j < gncell[1]; j++) {
	  idx[1] = j;
	  if (SqueezeDim[1]) idx[1] += NGH;
	  for (i = 0; i < gncell[0]; i++) {
	    idx[0] = i;
	    if (SqueezeDim[0]) idx[0] += NGH;
	    ms = idx[0]*strides[0]+idx[1]*strides[1]+idx[2]*strides[2];
	    m  = idx[0]*stride[0] +idx[1]*stride[1] +idx[2]*stride[2];
	    corr = -sv[l][m]-(p_int[ms]-p_int[ms+strides[l]])*		\
	      (inter[l][m]+inter[l][m+stride[l]])*invvol[m]/(2.*rho[m]);
	    sc[l]->field[ms] = corr;
	    /* Note that the test below is meant to avoid the
	       outermost ghost zone layer because we do not
	       necessarily know the source term there (in particular
	       because we cannot evaluate the potential derivative
	       there) */
	    if ((idx[l] > 0) && (idx[l] < gncell[l]-1) &&		\
		((gncell[lp1] == 1) || ((idx[lp1] > 1) && (idx[lp1] < gncell[lp1]-1))) && \
		((gncell[lp2] == 1) || ((idx[lp2] > 1) && (idx[lp2] < gncell[lp2]-1)))) {
	      corr = fabs(corr) / fabs (sv[l][m]);
	      if ((idx[l] > NGH) && (idx[l] < gncell[l]-NGH) &&		\
		  ((gncell[lp1] == 1) || ((idx[lp1] > NGH) && (idx[lp1] < gncell[lp1]-NGH))) && \
		  ((gncell[lp2] == 1) || ((idx[lp2] > NGH) && (idx[lp2] < gncell[lp2]-NGH)))) {
		rnbcorr++;
		rcorrmean += corr;
		if (rcorrmin > corr) rcorrmin = corr;
		if (rcorrmax < corr) rcorrmax = corr;
	      }
	      nbcorr++;
	      corrmean += corr;
	      if (corrmin > corr) corrmin = corr;
	      if (corrmax < corr) corrmax = corr;
	    }
	  }
	}
      }
      corrmean /= (real)nbcorr;
      rcorrmean /= (real)rnbcorr;
      pInfo ("\nStatistics of |source corr|/|source|, grid #%ld at level %ld,\n",\
	     fluid->desc->number, fluid->desc->level);
      pInfo ("for the correction in %s:\n", CoordNames[CoordType*3+InvCoordNb[l]]);
      pInfo ("Minimum : %lf\tActive mesh only %lf\n", corrmin, rcorrmin);
      pInfo ("Maximum : %lf\tActive mesh only %lf\n", corrmax, rcorrmax);
      pInfo ("Average : %lf\tActive mesh only %lf\n", corrmean, rcorrmean);
    }
  }
  CurrentToPatch (fluid); 
}
