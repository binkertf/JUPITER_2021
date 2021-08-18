/****************************************************/
/*                                                  */
/*                                                  */
/*  Fake MPI functions library for sequential built */
/*                                                  */
/*                                                  */
/****************************************************/

#include <stdio.h>
#include "mpi_dummy.h"
#include <unistd.h>
#include <sys/times.h>
void pInfo (const char *template, ...);

void MPI_Comm_rank (a, b)
     int a;
     int *b;
{
  *b = 0;			/* Only one process, with rank zero... */
}

void MPI_Comm_size (a, b)
     int a;
     int *b;
{
  *b = 1;			/* Only one process in the world communicator... */
}

void MPI_Init (argc, argv)
     int *argc;
     char **argv[];
{
  int test;
  pInfo ("\n       !!!! NOTE !!!!\n\n");
  pInfo ("This is a sequential built of the %s code\n", *argv[0]);
  pInfo ("If you planned to run the MPI-parallel version,\n");
  pInfo ("then you MUST rebuild the executable. Go into the\n");
  pInfo ("source directory (normally src/), then issue:\n");
  pInfo ("\ngmake BUILD=parallel\n");
  pInfo ("\nAny further invocation of gmake will refer to a parallel built.\n");
  MPI_Comm_size (0, &test);
  if (test != 1) {
    fprintf (stderr, "\n\nFurther warning : this sequential version has been compiled with a\n");
    fprintf (stderr, "wrong return value of the fake MPI_Comm_size : %d instead of 1\n",test);
    fprintf (stderr, "This trick is used in debugging or development phase only !!!\n");
    fprintf (stderr, "If you find strange results edit src/mpi_dummy.c and recompile\n");
  }
}

void MPI_Finalize ()
{
  int test;
  MPI_Comm_size (0, &test);
  if (test != 1) {
    fprintf (stderr, "\nWarning : this sequential version has been compiled with a\n");
    fprintf (stderr, "wrong return value of the fake MPI_Comm_size : %d instead of 1\n",test);
    fprintf (stderr, "If you find strange results edit src/mpi_dummy.c and recompile\n");
  }
}

void MPI_Bcast ()
{
}

void MPI_Isend ()
{
}

void MPI_Irecv ()
{
}

void MPI_Send ()
{
}

void MPI_Recv ()
{
}

void MPI_Barrier ()
{
}

void MPI_Wait ()
{
}

void MPI_Waitall ()
{
}

void MPI_Allreduce (void *ptr, void *ptr2, int count, int type, int foo3, int foo4)
{
  int i;
  for (i = 0; i < count; i++) {
    switch (type) {
    case MPI_DOUBLE:
      *(((double *)ptr2)+i) = (double)(*(((double *)ptr)+i));
      break;
    case MPI_INT:
      *(((int *)ptr2)+i) = (int)(*(((int *)ptr)+i));
      break;
    }
  }
}

double MPI_Wtime ()
{
  struct tms buffer;
  long Ticks;
  Ticks = sysconf (_SC_CLK_TCK);
  times (&buffer);
  return ((double)buffer.tms_utime/(double)Ticks);

}
