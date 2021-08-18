#include "jupiter.h"

void GridBuild (grids)
     GridFileInfo *grids;
{
  char *numbname[]={"first","second","third"};
  long i,level,levmax;
  long j,k,l,cellsize;
  tGrid *grid;
  real xmin[3], xmax[3], frac;
  i=0;
  do {
    grid = prs_malloc (sizeof(tGrid));
    level = grid->level = grids[i].level;
    grid->number = grids[i].number;
    grid->linenumber = grids[i].linenumber;
    for (j= 0; j < 3; j++) {
      xmin[j] = grids[i].xmin[j];
      xmax[j] = grids[i].xmax[j];
      grid->corner_min[j] = xmin[j];
      grid->corner_max[j] = xmax[j];
      grid->monoCPU = grids[i].monoCPU;
      grid->size[j] = grids[i].xmax[j]-grids[i].xmin[j];
      if ((grids[i].xmax[j] > corner_max0[j]) && (Periodic[j])) {
	pWarning ("Mesh %ld is beyond the upper template\n", grid->number); 
	pWarning ("limit along periodic dim %ld.\n", i);
	pWarning ("I assume that you know what you are doing...\n");
      }
      if ((grids[i].xmin[j] < corner_min0[j]) && (Periodic[j])) {
	pWarning ("Mesh %ld is beyond the lower template\n", grid->number); 
	pWarning ("limit along periodic dim %ld.\n", i);
	pWarning ("I assume that you know what you are doing...\n");
      }
      if ((grids[i].xmin[j] > grids[i].xmax[j])) {
	pWarning ("Mesh %ld has a lower limit (%g)", grid->number, grids[i].xmin[j]);
	pWarning ("larger than the upper limit (%g)!", grids[i].xmax[j]);
	pWarning ("(for dimension %ld)", j);
	prs_stderr ("Ill-ordered refinement information ? (see %s/warning.log)\n", OUTPUTDIR);
	prs_exit (1);
      }
      if (((grids[i].xmin[j] == grids[i].xmax[j])) && (j < NDIM)) {
	pWarning ("Mesh %ld has a lower limit (%g)", grid->number, grids[i].xmin[j]);
	pWarning ("identical to the upper limit (%g)!", grids[i].xmax[j]);
	pWarning ("(for dimension %ld)", j);
      }
      if ((grid->size[j] > (corner_max0[j]-corner_min0[j])) && (Periodic[j])) {
	prs_stderr ("Error: mesh %ld has size %g along dimension %ld,\n",\
		    grid->number, grid->size[j], i);
	prs_error  ("which is periodic with period: %g (smaller)\n", \
		    corner_max0[j]-corner_min0[j]);
      }
      for (k = INF; k <= SUP; k++) {
	if (Periodic[j] && grids[i].bc[j+k*3] && (j < NDIM)) {
	  pWarning ("You have specified a true boundary");
	  pWarning ("condition on a periodic dimension");
	  pWarning ("for grid at line %ld, on %s dimension, %s face",\
		      grids[i].linenumber, numbname[j], (k ? "upper": "lower"));
	  pWarning ("You have specified the BC code %ld.", grids[i].bc[j+k*3]);
	  pWarning ("I assume you meant 0.");
	  grids[i].bc[j+k*3] = 0;
	}
	grid->BoundaryConditions[j][k] = grids[i].bc[j+k*3];
      }
      l = level;
      levmax = LevMax;
      if (!Refine[j]) { 
	levmax = l = 0; 
      }
      grid->ncorner_min[j] = grids[i].nc_min[j];
      grid->ncorner_max[j] = grids[i].nc_max[j];
      grid->nsize[j] = grid->ncorner_max[j]-grid->ncorner_min[j];
      cellsize = 1 << (levmax-l+1);
      grid->dn[j] = cellsize;
      grid->ncell[j] = grid->nsize[j]/cellsize;
      grid->gncell[j] = grid->ncell[j]+2*Nghost[j];
      grid->gncorner_min[j] = grid->ncorner_min[j]-Nghost[j]*(1 << (levmax+1-l));
      grid->gncorner_max[j] = grid->ncorner_max[j]+Nghost[j]*(1 << (levmax+1-l));
      grid->gnsize[j] = grid->gncorner_max[j]-grid->gncorner_min[j];
    }
    for (j = 0; j < 3; j++) {
      grid->Edges[j] = (real *)prs_malloc((grid->gncell[j]+1)*sizeof(real));
      for (k = 0; k <= grid->gncell[j]; k++) {
	frac = (real)(k*grid->dn[j]+grid->gncorner_min[j])/\
	  (real)(ncorner_max0[j]-ncorner_min0[j]);
	grid->Edges[j][k] = frac*(corner_max0[j]-corner_min0[j])+corner_min0[j];
      }
    }
    grid->minCPU = 0;
    grid->maxCPU = CPU_Number-1; /* temporary cpu repartition */
    grid->next = GridList; /* We insert the newly created grid in the chained list of grids */
    grid->prev = NULL;

    if (GridList != NULL)
      GridList->prev = grid;
    GridList = grid;
  } while (!grids[i++].last);
}  

    
