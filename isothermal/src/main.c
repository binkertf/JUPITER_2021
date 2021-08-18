#include "jupiter.h"

int main (argc,argv)
     int argc;
     char *argv[];
{
  long  NbRestart, i, TimeStepCount=0, Iteration=0;
  real NextDate;
	real val;
  mpiwarning ();
  MPI_Init (&argc, &argv);
  MPI_Comm_rank (MPI_COMM_WORLD, &CPU_Rank);
  MPI_Comm_size (MPI_COMM_WORLD, &CPU_Number);
  NbRestart = GlobalInit (argc, argv);
  ResetChronometer (2);
  if (SmoothTaper == YES) GlobalDateInit=GlobalDate;
//  printf("GlobalDateInit in main: %lg\n",GlobalDateInit);
    for (i = NbRestart*NINTERM; i <= NTOT; i++) {
    ResetChronometer (1);
    if (MonitorCons) MonitorConservative ();
    if (TorqueMon && i>NbRestart*NINTERM ){
      MonitorTorque();
      if (NDIM == 3) MonitorTorqueZ (); 
    }
    if (KineticMon) MonitorKineticEnergy();
    NextDate = GlobalDate + DT;
    if (i % NINTERM == 0) {
      Write (CurrentOutputNumber++);
      if (Disable) prs_error("Run disabled.");
    }
    Iteration = 0;
    while (GlobalDate < NextDate-GlobalDate/1.e12) {
      pInfo ("[@]%g\t%ld\t", GlobalDate, CurrentOutputNumber);
      pInfo ("%ld\t%ld\t", i, TimeStepCount++);
      if ((!QuietOutput) && (!CPU_Rank)) {
	printf ("[%5.1f%%] of DT #%ld. Iter# %ld. ",
		(1.-(NextDate-GlobalDate)/DT)*100.0, i+1, ++Iteration);
	if ((i/NINTERM+1)*NINTERM <= NTOT)
	  printf ("Output #%ld after DT #%ld\r",\
		  CurrentOutputNumber, (i/NINTERM+1)*NINTERM);
	else
	  printf ("No other output scheduled\r");
	fflush (stdout);
      }
      ResetChronometer (0);
//	printf(" dates %.17g %.17g %.17g %.17g\n",GlobalDate,NextDate,NextDate-1e-12,NextDate-GlobalDate); 
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
