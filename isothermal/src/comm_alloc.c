#include "jupiter.h"

static long CommNb=0;
MPI_Request *Req=NULL;

void Comm_Alloc (start)
     Communicator *start;
{
  Communicator *com;
  long size,i,nvar;
  long cpusrc, cpudest;
  com = start;
  while (com != NULL) {
    cpusrc  = com->CPU_src; /* Beware ! size of com is always size of dest... */
    cpudest = com->CPU_dest;
    if ((cpusrc == CPU_Rank) || (cpudest == CPU_Rank)) {
      if (cpusrc != cpudest)
	CommNb++;
      size = 1;
      for (i = 0; i < NDIM; i++)
	size *= (com->imax_dest[i]-com->imin_dest[i]);
      com->size = size;
      nvar = 1+NDIM;		/* Density (or mass flux) + velocities (or momentum flux) */
      nvar++;	/* Energy or energy flux */
      if (com->type == FLUX)
	    nvar++;			/* pressure information */

  //if (DUSTDIFF==TRUE)
    //nvar++;     //diffusion flux

      /* It needs to be allocated even in the isothermal case for the first communication */
      com->buffer = prs_malloc (nvar*size*sizeof(real)*NbFluids);
    }
    com = com->next;
  }
  Req = (MPI_Request *)realloc(Req, (size_t)((CommNb+1)*sizeof(MPI_Request)));
  /* The "+1" is here to prevent Electric Fence considering the above instruction as */
  /* a bug when running only on one process */
}

long dimfield (field)
     long field;
{
  long d;
  switch (field) {
  case _Velocity_:
    d=NDIM;
    break;
  default:
    d = 1;
  }
  return d;
}

long nbvariables (size, fieldtype)
     long size;
     int *fieldtype;
{
  long i, nvar = 0;
  for (i = 0; i < size; i++) {
    nvar += dimfield(fieldtype[i]);
  }
  return nvar;
}
