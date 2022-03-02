#include "jupiter.h"

int main (argc,argv)
     int argc;
     char *argv[];
{
  long  NbRestart, i, TimeStepCount=0, Iteration=0;
  real NextDate;
  mpiwarning ();
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  NbRestart = GlobalInit (argc, argv);
  ResetChronometer (2);
  if (SmoothTaper == YES) GlobalDateInit=GlobalDate;
  for (i = NbRestart*NINTERM; i <= NTOT; i++) {
    ResetChronometer (1);
    if (MonitorCons) MonitorConservative ();
    if (TorqueMon && i>NbRestart*NINTERM ){
      MonitorTorque();
      if (NDIM == 3) MonitorTorqueZ ();
    }
    //if (EneMon) 
    MonitorEnergy();
    if (DUSTSUBL){
      if (i==NbRestart*NINTERM) ReadSublDustMass(NbRestart);
      MonitorSublimatedDustMass ();
    }
    if ((!DUSTUCAP)==NO){
      if (i==NbRestart*NINTERM) ReadRemDustMass(NbRestart);
      MonitorRemovedDustMass ();
    } 
    NextDate = GlobalDate + DT;
    if (i % NINTERM == 0) {
//      CommAll ();
      Write (CurrentOutputNumber++);
      if (Disable) prs_error("Run disabled.");
    }
    Iteration = 0;
    while (GlobalDate < NextDate-GlobalDate/1.e12) {
      pInfo ("[@]%g\t%ld\t", GlobalDate, CurrentOutputNumber);
      pInfo ("%ld\t%ld\t", i, TimeStepCount++);
      if (!QuietOutput) {
	if (!CPU_Rank)
	  printf ("[%5.1f%%] of DT #%ld. Iter# %ld. ",
		  (1.-(NextDate-GlobalDate)/DT)*100.0, i+1, ++Iteration);
	if ((i/NINTERM+1)*NINTERM <= NTOT) {
	  if (!CPU_Rank)
	    printf ("Output #%ld after DT #%ld\r",		\
		    CurrentOutputNumber, (i/NINTERM+1)*NINTERM);
	} else {
	  MPI_Finalize ();
	    if (!CPU_Rank)
	      printf ("End of run !\n");
	    return 0;
	}
	fflush (stdout);
      } else {
	printf (".");
	fflush (stdout);
      }
      ResetChronometer (0);
      GlobalDate += RecursiveIteration (NextDate-GlobalDate, 0L);
      ReadChronometer (0, "one full hydro time step");
    }
    prs_msg ("[100.0%%] of DT #%ld.                                          \n", i+1);
    fflush (stdout);
    ReadChronometer (1, "one full DT integration");
  }
  ReadChronometer (2, "full run (initialization phase excepted)");
  MPI_Finalize ();
  return 0;
}
