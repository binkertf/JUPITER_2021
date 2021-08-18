#include "jupiter.h"
extern MPI_Request *Req;

void ExecCommOneVar (long, long, real *); 
void ExecCommSameOneField (long, real *); 


void ExecCommSameOneField (long lev, real *myfield) 
{	
  /* Execute intra-level communications, for level lev. These are
     GHOST type communications */
  Communicator *com;
  long i,j,k,l,m,n,imin[3],imax[3],stride[3];
  FluidPatch *fluid;
  real *source[80];
  com = ComListGhost;
  while (com != NULL) {
    if ((com->dest_level == lev) && (com->src_level == lev)) {
      if (com->CPU_src == CPU_Rank) {
	l=0; 
	fluid = com->srcg->Fluid;
	while (fluid != NULL) {
	  source[l++] = myfield;
		fluid = fluid->next;
	}
	for (l = 0; l < 3; l++) {
	  imin[l] = com->imin_src[l];
	  imax[l] = com->imax_src[l];
	  stride[l] = com->srcg->stride[l];
	}
	n=0;
	for (l=0; l < NbFluids; l++) { 					
	  for (k = imin[2]; k < imax[2]; k++) {
	    for (j = imin[1]; j < imax[1]; j++) {
	      for (i = imin[0]; i < imax[0]; i++) {
		m = i*stride[0]+j*stride[1]+k*stride[2];
		com->buffer[n++] = source[l][m];
	      }
	    }
	  }
	}
      }
    }
    com = com->next;
  }
  ExecCommOneVar (lev,lev,myfield);
}


void ExecCommOneVar (long levsrc, long levdest, real *myfield) 
{
  Communicator *com=NULL, *start=NULL;
  long i,nbreq=0,l,imin[3],imax[3],stride[3];
  long j,k,m,n,sqz[3],ii[3];
  real *dest[80];
  FluidPatch *fluid;
  MPI_Status stat;
  tGrid_CPU *desc;
  start = com = ComListGhost;
  while (com != NULL) {		/* We initiate non-blocking communications */
    if ((com->dest_level == levdest) &&\
	(com->src_level == levsrc) &&\
	(com->type == GHOST)) {
      if (com->CPU_src != com->CPU_dest) {
       	if (com->CPU_src == CPU_Rank) {
	  MPI_Isend (com->buffer, com->size*NbFluids, MPI_DOUBLE,\
		     com->CPU_dest, com->tag, MPI_COMM_WORLD, Req+nbreq);
	  nbreq++;
	}
	if (com->CPU_dest == CPU_Rank) {
	  MPI_Irecv (com->buffer, com->size*NbFluids, MPI_DOUBLE,\
		     com->CPU_src, com->tag, MPI_COMM_WORLD, Req+nbreq);
	  nbreq++;
	}
      }
    }
    com=com->next;
  }
  for (i = 0; i < nbreq; i++)
    MPI_Wait (Req+i, &stat);
  com = start;
  while (com != NULL) {		/* We use the buffers to update the
				   ghost or "mean" values */
    if ((com->dest_level == levdest) &&		\
	(com->src_level == levsrc) &&		\
	(com->type == GHOST)) {
      if (com->CPU_dest == CPU_Rank) {
	l=0; 
	fluid = com->destg->Fluid;
	while (fluid != NULL) {
	  dest[l++] = myfield;
	  fluid = fluid->next;
	}
	desc = com->destg;
	for (l = 0; l < 3; l++) {
	  imin[l] = com->imin_dest[l];
	  imax[l] = com->imax_dest[l];
	  stride[l] = com->destg->stride[l];
	  if (HydroStaticReady)
	    sqz[l] = desc->Fluid->Rho_eq_c->stride[l];
	}
	n=0;
	for (l=0; l < NbFluids; l++) { 											/* here we run over all the fluids */
	  for (k = imin[2]; k < imax[2]; k++) {
	    for (j = imin[1]; j < imax[1]; j++) {
	      for (i = imin[0]; i < imax[0]; i++) {
		ii[0] = i;
		ii[1] = j;
		ii[2] = k;
		m = i*stride[0]+j*stride[1]+k*stride[2];
		dest[l][m] = com->buffer[n++];
	      }
	    }
	  }
	}
      }
    }
    com = com->next;
  }
}
