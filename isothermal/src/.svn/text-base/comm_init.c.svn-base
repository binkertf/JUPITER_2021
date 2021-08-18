#include "jupiter.h"

extern long cpugrid_number;

struct destblocks {
  long block[27][6];
  long shift[27][3];
  long nb;
};

typedef struct destblocks DestBlocks;

/* This file implements the communicators list build up algorithm */


/* The following function splits a block into pieces if it happens to
have parts that lie outside of the main domain along each dimension
for which the system is periodic. It returns as many pieces (27 at
maximum) that all lie within the main domain, and that have a
one-to-one correspondance with the constitutive blocks of the initial
cube, modulo a translation by an integer number of periods along the
relevant dimensions. */
inline void BlockButcher (cube, dblocks)
     long *cube;
     DestBlocks *dblocks;
{
  long i,j,k;
  boolean cut; 
  for (i = 0; i < 6; i++) {
    dblocks->block[0][i] = cube[i];
    dblocks->shift[0][i] = 0;
  }
  dblocks->nb = 1;		/* Fast track for most blocks */
  if ((cube[0] >= 0) && (cube[3] <= ncorner_max0[0]) &&\
      (cube[1] >= 0) && (cube[4] <= ncorner_max0[1]) &&	\
      (cube[2] >= 0) && (cube[5] <= ncorner_max0[2]))
    return;			/* Most cubes returns here. They do not overlap a periodic BC boundary. */
  do {
    cut = FALSE;
    for (i = 0; i < dblocks->nb; i++) {
      for (j = 0; j < NDIM; j++) {
	if (Periodic[j]) {
	  if ((dblocks->block[i][j] < 0) && (dblocks->block[i][j+3] > 0)) { /* This block sits on a lower boundary */
	    cut = TRUE;		/* it therefore needs to be cut in two parts */
	    for (k = 0; k < 3; k++) { /* We clone it */
	      dblocks->block[dblocks->nb][k] = dblocks->block[i][k];
	      dblocks->block[dblocks->nb][k+3] = dblocks->block[i][k+3];
	      dblocks->shift[dblocks->nb][k] = dblocks->shift[i][k];
	    }
	    dblocks->block[i][j] = 0; /* We modify the original. It lies within the main domain */
	    dblocks->block[dblocks->nb][j+3] = 0; /* The copy, that has rank nb+1, corresponds to the other bit */
	    dblocks->block[dblocks->nb][j] += ncorner_max0[j]; /* We shift it so that it lies within the main domain */
	    dblocks->block[dblocks->nb][j+3] += ncorner_max0[j];
	    dblocks->shift[dblocks->nb][j] += ncorner_max0[j]; /* We update the shift vector to remember where the original was. */
	    dblocks->nb++;
	  } else if ((dblocks->block[i][j] < ncorner_max0[j]) && (dblocks->block[i][j+3] > ncorner_max0[j])) {
	    cut = TRUE;		/* This block sits on an upper boundary */
	    for (k = 0; k < 6; k++) { /* it therefore needs to be cut in two parts */
	      dblocks->block[dblocks->nb][k] = dblocks->block[i][k]; /* We do similar things as before */
	      dblocks->block[dblocks->nb][k+3] = dblocks->block[i][k+3];
	      dblocks->shift[dblocks->nb][k] = dblocks->shift[i][k];
	    }
	    dblocks->block[i][j+3] = ncorner_max0[j]; /* The original lies within the main grid */
	    dblocks->block[dblocks->nb][j] = ncorner_max0[j]; /* The copy lies outside */
	    dblocks->block[dblocks->nb][j] -= ncorner_max0[j];
	    dblocks->block[dblocks->nb][j+3] -= ncorner_max0[j]; /* It is translated to have it within the main domain */
	    dblocks->shift[dblocks->nb][j] += -ncorner_max0[j];	/* We keep track of the shift... */
	    dblocks->nb++;
	  } else if (dblocks->block[i][j] < 0) { /* the block does not need splitting but it needs to be shifted */
	    cut = TRUE;		/* Note: we set here cut to TRUE (and
				   in the next test as well), although
				   nothing has been cut, simply
				   because the block thus shifted may
				   still lie outside of the main
				   domain. This happens only in very
				   special cases (not very common and
				   probably not even useful), e.g. if
				   one has a size of one along the
				   corresponding dimension, then if
				   NGH>1 some remote ghost zones may
				   well never be sent to the main
				   domain. */
	    dblocks->block[i][j] += ncorner_max0[j]; /* so we do... */
	    dblocks->block[i][j+3] += ncorner_max0[j];
	    dblocks->shift[i][j] += ncorner_max0[j]; /* ...and we keep track of the shift */
	  } else if (dblocks->block[i][j+3] > ncorner_max0[j]) { /* the block does not need splitting but it needs to be shifted */
	    cut = TRUE;
	    dblocks->block[i][j] -= ncorner_max0[j]; /* Same thing as before but on the other side */
	    dblocks->block[i][j+3] -= ncorner_max0[j];
	    dblocks->shift[i][j] += -ncorner_max0[j];
	  }
	  if (dblocks->nb > 27) 
	    prs_error ("Internal error. Block split in more than 27 parts. Check code in comm_init.c");
	}
      }
    }
  }
  while (cut);			/* As long as we have undertaken any
				   action that implies cutting a
				   block, we recheck entirely the
				   list. */
}


/* The following function checks whether two cubes intersect, in which
   case it returns their intersection. A properly dimensionned array
   (inter) has always to be provided as an argument in case the
   intersection is non-empty */
inline boolean CubeIntersect (cube1, cube2, inter)
     long *cube1, *cube2, *inter;
{
  long k;
  if ((cube1[0] >= cube2[3]) || (cube2[0] >= cube1[3]) ||\
      (cube1[1] >= cube2[4]) || (cube2[1] >= cube1[4]) ||\
      (cube1[2] >= cube2[5]) || (cube2[2] >= cube1[5]))
    return FALSE;
      /* Note that the above test assumes that cube[i+3]>cube[i],
	 i.e. that coordinates are properly ordered. CubeIntersect
	 gives an erroneous result if this is not the case */
  for (k = 0; k < 3; k++) {
    inter[k]   = (cube1[k]   < cube2[k]   ? cube2[k]   : cube1[k]);
    inter[k+3] = (cube2[k+3] > cube1[k+3] ? cube1[k+3] : cube2[k+3]);
  }
  return TRUE;
}
    
/* The following function checks whether two cubes intersect, but it
   does not return anything apart from a boolean result, even if the
   intersection is non void. It actually checks whether the cubes
   intersect with each other or with any of their replicants, in each
   dimension for which the mesh is periodic (it actually keeps one
   cube fixed, and moves the other one back and forth one period in
   each periodic dimension. It returns as soon as it can give a
   positive answer (the cubes intersect). */
inline boolean CubeIntersectPeriodic (cube1, cube2)
     long *cube1, *cube2;
{
  long k[3], cube[6], l;
  long periodic[3];
  for (l = 0; l < 3; l++) {
    periodic[l] = Periodic[l];
    if ((periodic[l]) && (cube1[l] >= 0)  && (cube1[l+3] <= ncorner_max0[l]) &&\
	(cube2[l] >= 0)  && (cube2[l+3] <= ncorner_max0[l]))
      periodic[l] = 0;
  }
  /* The above test means that there is no need to consider replicants
     along a periodic dimension if both cubes lie within the main
     period. It significantly speeds up the communicator build up time when
     several dimensions are periodic */
  for (k[0] = -(periodic[0] == TRUE); k[0] < 1+(periodic[0] == TRUE); k[0]++) {
    for (k[1] = -(periodic[1] == TRUE); k[1] < 1+(periodic[1] == TRUE); k[1]++) {
      for (k[2] = -(periodic[2] == TRUE); k[2] < 1+(periodic[2] == TRUE); k[2]++) {
	for (l = 0; l < 6; l++)
	  cube[l] = cube1[l]+k[l%3]*ncorner_max0[l%3];
	if ((cube[0] < cube2[3]) && (cube2[0] < cube[3]) &&	\
	    (cube[1] < cube2[4]) && (cube2[1] < cube[4]) &&	\
	    (cube[2] < cube2[5]) && (cube2[2] < cube[5]))
	  return TRUE;
      }
    }
  }
  /* Note that the above test assumes that cube[i+3]>cube[i],
     i.e. that coordinates are properly ordered. CubeIntersectPeriodic
     gives an erroneous result if this is not the case */
  return FALSE;
}
    
void BuildCommunicators () {
  int TagNumber = 0;
  tGrid_CPU *destg, *srcg;
  Communicator *com;
  CommHash *hash;
  DestBlocks dblocks;
  long cdest[6], csrc[6], ccom[6];
  long i,j,dim,k,stripped, fImin[3], fImax[3];
  destg = Grid_CPU_list;
  pInfo ("Building communicators...\n");
  DestroyCom (&ComListGhost);		/* Destroy any previous list of communicators, if any */
  DestroyCom (&ComListFlux);
  DestroyCom (&ComListMean);
  /* Allocate hash table of chained lists of communicators that share
     the same source CPU grid, and initialize */
  CommHashSrc = (CommHash **)prs_malloc(sizeof(CommHash *)*(cpugrid_number));
  for (i = 0; i < cpugrid_number; i++)
    CommHashSrc[i] = NULL;
  while (destg != NULL) {
    srcg = Grid_CPU_list;
    while (srcg != NULL) {	/* The following test corresponds to
				   the creation of communicators that
				   will fill boundaries of level l,
				   coming from grids at level l
				   and l-1. It will also serve to
				   create flux information
				   communicators from level l to l-1.
				   Note that communicators between
				   level l-1 and l are considered only
				   if the corresponding iface value of
				   level l is -1, otherwise this means that
				   we are considering an inter-CPU face on a
				   mesh split by MPI. */
      if ((srcg->level == destg->level) || (srcg->level == destg->level-1)) {
	for (i = 0; i < 3; i++) {
	  cdest[i]   = destg->gncorner_min[i];
	  cdest[i+3] = destg->gncorner_max[i];
	  csrc[i]    = srcg->ncorner_min[i];
	  csrc[i+3]  = srcg->ncorner_max[i];
	}
	/* We check whether an intersection is possible */
	
	if (CubeIntersectPeriodic (cdest, csrc)) { /* We need to check
						      cubes+replicants
						      intersections */
	  for (dim = 0; dim < NDIM; dim++) {
	    for (j = INF; j <= SUP; j++) {
	      if ((srcg->level == destg->level) || (destg->iface[dim][j]==-1)) { 
		/* This test means that either we consider a communicator at the same
		   level (l --> l), or we consider a communicator from level l-1-->l,
		   but only on a face of cpugrid "dest" that corresponds to a face of
		   the corresponding parent grid, rather than an internal interprocess
		   interface, in which case there is neither a need to first fill the
		   ghosts with lower level interpolation, nor there is a need to send
		   the flux of conservative quantities to the lower level at the end
		   of a time step sequence on that level*/
		for (i = 0; i < 3; i++) {
		  cdest[i]   = destg->gncorner_min[i];
		  cdest[i+3] = destg->gncorner_max[i];
		}
		cdest[dim]   = (j == INF ? destg->gncorner_min[dim] : destg->ncorner_max[dim]);
		cdest[dim+3] = (j == INF ? destg->ncorner_min[dim] : destg->gncorner_max[dim]);
		BlockButcher (cdest, &dblocks);
		for (i = 0; i < dblocks.nb; i++) {
		  if (CubeIntersect (dblocks.block[i], csrc, ccom)) {
		    com = (Communicator *)prs_malloc (sizeof(Communicator));
		    for (k = 0; k < 3; k++) {
		      com->Imin[k]       = ccom[k]-dblocks.shift[i][k];	/* In the absolute coordinates we keep a track of the */
		      com->Imax[k]       = ccom[k+3]-dblocks.shift[i][k]; /* ORIGINAL position (may lie outside the main domain) */
		      /* The reason for that is that otherwise we
			 would strip far too many communicators in
			 case of several periodic dimensions, and we
			 would loose the corners or edges
			 communicators -they would be detected as
			 being included respectively in edges or faces
			 communicators- */
		      com->imin_src[k]   = (ccom[k]   -  srcg->gncorner_min[k])/\
			(srcg->Parent->dn[k]); 
		      com->imax_src[k]   = (ccom[k+3] -  srcg->gncorner_min[k])/\
			(srcg->Parent->dn[k]); 
		      com->imin_dest[k]  = (ccom[k]-dblocks.shift[i][k]   -  destg->gncorner_min[k])/\
			(destg->Parent->dn[k]); 
		      com->imax_dest[k]  = (ccom[k+3]-dblocks.shift[i][k] -  destg->gncorner_min[k])/\
			(destg->Parent->dn[k]); 
		    }
		    com->CPU_dest    = destg->cpu;
		    com->CPU_src     = srcg->cpu;
		    com->facedim     = dim;
		    com->faceside    = j;
		    com->grid_dest   = destg->parent;
		    com->grid_src    = srcg->parent;
		    com->srcg        = srcg;
		    com->destg       = destg;
		    com->dest_level  = destg->level;
		    com->src_level   = srcg->level;
		    com->nb_src      = srcg->number;
		    com->nb_dest     = destg->number;
		    com->tag         = TagNumber++;
		    com->type        = GHOST;
		    com->prev        = NULL;
		    com->next        = ComListGhost;
		    if (ComListGhost != NULL) ComListGhost->prev = com;
		    ComListGhost          = com;
		    hash = (CommHash *)prs_malloc(sizeof(CommHash)); /* We update the hash table */
		    hash->com = com;
		    hash->next = CommHashSrc[com->nb_src];
		    hash->prev = NULL;
		    if (hash->next != NULL) hash->next->prev = hash;
		    CommHashSrc[com->nb_src] = hash;
		    if (destg->level == srcg->level+1) {
		      /* We may then also create a flux communicator */
		      for (k = 0; k < 3; k++) {
			fImin[k] = ccom[k]-dblocks.shift[i][k];	/* fImin & fImax therefore stand for the communicator */
			fImax[k] = ccom[k+3]-dblocks.shift[i][k]; /* boundaries of the ghost which may lie outside of the */
			/* main part of the mesh, in the dimension(s) for which the BC are periodic. */
				/* This allows the direct comparison
				   below (destg stands for the ghost
				   com dest, ie the highest level of
				   the two) */
			if (k != dim) {	/* This is for dimensions perpendicular to the current direction */
			  if (destg->ncorner_min[k] > fImin[k]) fImin[k] = destg->ncorner_min[k];
			  if (destg->ncorner_max[k] < fImax[k]) fImax[k] = destg->ncorner_max[k];
			} else {
			  if (j == SUP)
			    fImax[k] = fImin[k]+1; /* fImin is at the upper face of the CPU grid */
			  else
			    fImin[k] = fImax[k]-1; /* fImax is at the lower face of the CPU grid */
			}
		      }
		      /* Would that communicator be non void ? Some
			 ghost communicators from level l-1=>l would
			 have no associated flux communicators from
			 level l=>l-1, such as edges and corners: */
/* Note the test on the second line below. It is here to avoid to
   create 2 flux communicators in case a ghost com was sitting on top
   of the (periodic) template edge. Only the block next to the upper
   level cpugrid is then taken into account. */
		      if ((fImin[0] < fImax[0]) && (fImin[1] < fImax[1]) && (fImin[2] < fImax[2])\
			  && ((fImin[dim] == destg->ncorner_max[dim]) || (fImax[dim] == destg->ncorner_min[dim]))) {
			com  = (Communicator *)prs_malloc (sizeof(Communicator)); /* OK it is non-void */
				/* Note that now the source is what
				   was the destination of the
				   associated ghost communicators,
				   while the destination was the
				   source. Hence the swap src <==>
				   dest in the affectations below. */
			for (k = 0; k < 3; k++) {
			  /* By construction the fImin & fImax are
			     relative to the (ghost) dest position,
			     hence the evaluation below: */
			  com->imin_src[k] =(fImin[k]-destg->gncorner_min[k])/(destg->Parent->dn[k]); 
			  com->imax_src[k] =(fImax[k]-destg->gncorner_min[k])/(destg->Parent->dn[k]); 
			  /* On the other hand the (ghost) src is
			     necessarily on the main domain, so the
			     communicator coordinates are shifted onto
			     this latter. This has the consequence
			     that the fluxes on the edge of a periodic
			     mesh will be corrected on the adequate
			     side (i.e. where the source of the ghost
			     com is)*/
			  com->imin_dest[k]=(fImin[k]+dblocks.shift[i][k]-srcg->gncorner_min[k])/(srcg->Parent->dn[k]); 
			  com->imax_dest[k]=(fImax[k]+dblocks.shift[i][k]-srcg->gncorner_min[k])/(srcg->Parent->dn[k]);
			  if (k == dim) { 
			    /* The expression below works either on
			       INF or SUP faces of a CPU grid. It
			       obviously works on a SUP side, while on
			       an INF side the integer division above
			       (for imin) was "floored" to
			       "edge"-1  */
			    com->imax_dest[k] = com->imin_dest[k]+1;
			    com->imax_src[k] = com->imin_src[k]+1;
			  }
			  com->Imin[k] = fImin[k];
			  com->Imax[k] = fImax[k];
			}
			com->CPU_dest    = srcg->cpu; /* Source and destination are swapped wrt ghost communicator */
			com->CPU_src     = destg->cpu;
			com->facedim     = dim;
			com->faceside    = j;
			com->grid_dest   = srcg->parent; /* same thing */
			com->grid_src    = destg->parent;
			com->dest_level  = srcg->level; /* same thing */
			com->src_level   = destg->level;
			com->srcg        = destg;
			com->destg       = srcg;
			com->nb_src      = destg->number;	/* same thing */
			com->nb_dest     = srcg->number;
			com->tag         = TagNumber++;
			com->type        = FLUX;
			com->prev        = NULL;
			com->next        = ComListFlux;
			if (ComListFlux != NULL) ComListFlux->prev = com;
			ComListFlux      = com;
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      if (srcg->level == destg->level+1) { /* We look for MEAN communicators */
	for (i = 0; i < 3; i++) {
	  cdest[i]   = destg->ncorner_min[i];
	  cdest[i+3] = destg->ncorner_max[i];
	  csrc[i]    = srcg->ncorner_min[i];
	  csrc[i+3]  = srcg->ncorner_max[i];
	}
	if (CubeIntersectPeriodic (cdest, csrc)) { 
	  /* This test is fine there, since both cdest and csrc are
	     ghost excluding, so they lie on the main domain. However
	     one could write a faster test */
	  for (dim = 0; dim < NDIM; dim++) {
	    for (j = INF; j <= SUP; j++) {
	      if (srcg->iface[dim][j]==-1) {
		for (i = 0; i < 3; i++) {
		  csrc[i]   = srcg->ncorner_min[i];
		  csrc[i+3] = srcg->ncorner_max[i];
		}
		csrc[dim]   = (j == INF ? srcg->ncorner_min[dim] : srcg->ncorner_max[dim]-destg->Parent->dn[dim]);
		csrc[dim+3] = (j == INF ? srcg->ncorner_min[dim]+destg->Parent->dn[dim] : srcg->ncorner_max[dim]);
		if (CubeIntersect (cdest,csrc,ccom)) {
		  com = (Communicator *)prs_malloc (sizeof(Communicator));
		  for (k = 0; k < 3; k++) {
		    com->Imin[k]      = ccom[k];
		    com->Imax[k]      = ccom[k+3];
		    com->imin_src[k]  = (ccom[k]   -srcg->gncorner_min[k])/(srcg->Parent->dn[k]); 
		    com->imax_src[k]  = (ccom[k+3] -srcg->gncorner_min[k])/(srcg->Parent->dn[k]); 
		    com->imin_dest[k] = (ccom[k]   -destg->gncorner_min[k])/(destg->Parent->dn[k]); 
		    com->imax_dest[k] = (ccom[k+3] -destg->gncorner_min[k])/(destg->Parent->dn[k]); 
		  }
		  com->CPU_dest    = destg->cpu;
		  com->CPU_src     = srcg->cpu;
		  com->facedim     = dim;
		  com->faceside    = j;
		  com->grid_dest   = destg->parent;
		  com->grid_src    = srcg->parent;
		  com->srcg        = srcg;
		  com->destg       = destg;
		  com->dest_level  = destg->level;
		  com->src_level   = srcg->level;
		  com->nb_src      = srcg->number;
		  com->nb_dest     = destg->number;
		  com->tag         = TagNumber++;
		  com->type        = MEAN;
		  com->prev        = NULL;
		  com->next        = ComListMean;
		  if (ComListMean != NULL) ComListMean->prev = com;
		  ComListMean          = com;
		}
	      }
	    }
	  }
	}
      }
      srcg = srcg->next;
    }
    destg = destg->next;
  }
  pInfo ("%d communicators found\n", TagNumber);
  SetOverlapFlag ();
  SetCommSource ();
  stripped = StripCommunicators (); /* This must be done _after_ setting Src & Dest flags*/
  pInfo ("%d communicators left\n", TagNumber-stripped);
  DestroyHash() ;		/* Clean up hash table */
  Comm_Alloc (ComListGhost);
  Comm_Alloc (ComListFlux);
  Comm_Alloc (ComListMean);
  CheckBC ();
}

