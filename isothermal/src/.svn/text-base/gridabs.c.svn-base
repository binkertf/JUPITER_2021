#include "jupiter.h"

void GridAbs (grids)
     GridFileInfo *grids;
{
  long i,l;
  i=0;
  do {
    if (grids[i].level > LevMax)
      LevMax = grids[i].level;
  } while (!grids[i++].last);
  pInfo ("Maximum grid level is %ld\n", LevMax);
  for (i = 0; i < 3; i++) {
    l = LevMax+1;
    if (!Refine[i]) l=1;
    ncorner_min0[i] = 0;
    ncorner_max0[i] = Ncell0[i]*(1<<l);
  }
}

void GridPos (grids)
     GridFileInfo *grids;
{
  long i,j,levmax,lp,level;
  i=0;
  do {
    level = grids[i].level;
    for (j = 0; j < 3; j++) {
      lp = (level > 0 ? level : 1);
      levmax = LevMax;
      if (!Refine[j]) { 
	levmax = 0; 
	lp = 1;
      }
      if (j < NDIM) {
	grids[i].nc_min[j] = (long)(((real)Ncell0[j]*(grids[i].xmin[j]-corner_min0[j])/\
				     (corner_max0[j]-corner_min0[j]))*\
				    (1<<(lp-1)))*(1<<(levmax+2-lp));
	grids[i].nc_max[j] = (long)(((real)Ncell0[j]*(grids[i].xmax[j]-corner_min0[j])/\
				     (corner_max0[j]-corner_min0[j]))*\
				    (1<<(lp-1))+.99999)*(1<<(levmax+2-lp));
      } else {
	grids[i].nc_min[j] = 0L;
	grids[i].nc_max[j] = 2L;
      }
    }
  } while (!grids[i++].last);
}

