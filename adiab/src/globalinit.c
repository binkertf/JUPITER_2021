#include "jupiter.h"

long GlobalInit (argc, argv)
     int argc;
     char *argv[];
{
      prs_msg ( "\n\n");
      prs_msg ( "                 .        ___---___                    . \n");                  
      prs_msg ( "        .              .--\\       --.     .     .         .\n");
      prs_msg ( "                     ./.;_.\\     __/~ \\.     \n");
      prs_msg ( "                    /;  / `-'  __\\    . \\              \n");              
      prs_msg ( "  .        .       / ,--'     / .   .;   \\        |\n");
      prs_msg ( "                  | .|       /       __   |      -O-       .\n");
      prs_msg ( "                 |__/    __ |  . ;   \\ | . |      |\n");
      prs_msg ( "                 |      /  \\_    . ;| \\___|    \n");
      prs_msg ( "    .    o       |      \\  .~\\___,--'     |           .\n");
      prs_msg ( "                  |     | . ; ~~~~\\_    __|\n");
      prs_msg ( "     |             \\    \\   .  .  ; \\  /_/   .\n");
      prs_msg ( "    -O-        .    \\   /         . |  ~/                  .\n");
      prs_msg ( "     |    .          ~\\ \\   .      /  /~          o \n");
      prs_msg ( "   .                   ~--___ ; ___--~       \n");
      prs_msg ( "                  .          ---         .         \n");     
      prs_msg ( "             __  _  _  ____  __  ____  ____  ____ \n");  
      prs_msg ( "           _(  )/ )( \\(  _ \\(  )(_  _)(  __)(  _ \\ \n");  
      prs_msg ( "          / \\) \\) \\/ ( ) __/ )(   )(   ) _)  )   /\n");  
      prs_msg ( "          \\____/\\____/(__)  (__) (__) (____)(__\\_)\n");       
      prs_msg ( "\n\n");

  long NbRestart=0;
  char parameters[MAXLINELENGTH];
  GridFileInfo grids[MAXGRIDS];
  tGrid *grid;
  srand (CPU_Rank);
  ReadDefaultInOut (); 
  NbRestart = switches (argc, argv, parameters); /* We parse the command line */
  /* This gives us 'parameters', the name of the parameter file */
  ReadVariables (parameters);	/* which we parse. */
  SubsDef ("output directory", OUTPUTDIR, DefaultOut);
  FlushLog ();			/* Since we now know where OUTPUTDIR
				   is, we can write the log files
				   (info, warning and error), see
				   log.c */
  checkpara ();			/* We check whether the user has not
				   mistakingly started a sequential
				   build with an MPI command */
  DumpInitCode ("init.log"); /* We write the initial conditions to 'init.log' */
  findcodepot (POTENTIALCODE, &PotentialCode); /* We check whether we know the potential code */
  DumpPotCode ("pot.log"); /* We write the potential expression to 'pot.log' */
  DumpParameters ();		/* Write input parameters to output directory */
  DumpSources (argc, argv);	/* Write hidden tarball source file to
				   output directory */
  TestEndian ();		/* This is to let IDL knows whether to swap
				   endians during post-processing */
  if (Merge) merge(NbRestart);
  if (AddSubPatch) refine();	/* Everything is defined in SubPatchInfo */
  if (!Restart)
    ScanGridFile (GRIDFILE);
  else {
    ReadGrids (NbRestart, grids);
    ConstructGrids (grids);
    CurrentOutputNumber = NbRestart;
  }
  grid = GridList;
  while (grid != NULL) {splitgrid (grid); grid = grid->next;}
  BuildCommunicators ();
  if (KEPLERIAN && MERIDIANCORRECTION) {		/* Note that
							   calling the function below is
							   correct even in the case of a
							   restart because it relies upon the
							   analytical expression of the
							   initial conditions */
    ResetChronometer (3);
    GetMaxLevProjectedMeridian ();
    ForAllPatches (InitFluxCorrectionMeridian);
    ReadChronometer (3, "Full meridian flux correction initialization");
    /* Note that the above must be done EVEN IF MERIDIANCORRECTION is
       disabled, since it is used in the Hydrostatic Enforcement
       initialization */
  }
  if (HydroStaticEnforce) HydroStatPrepare ();
  InitWholeHierarchy (&NbRestart);
  if (HydroStaticEnforce) HydroStaticReady = YES;
  /* The above boolean regulates the choice of communication adaptation
     for ghost zone fillings (see comm_adapt.c). Note that since the
     initialization loops are performed only on the active grid, ghost
     zones are initialized afterwards via communications. In
     particular, the 'Energy' field is initialized via communications
     (which are set temporarily to non-isothermal so that this field
     is added to communication buffers) and this *must* be done the
     same way for the HSE initialization and the run
     initialization. This is why the flag is set to true only after
     both initialization have been performed */
  if (RayTracing) Ray_Tracing ();
  return NbRestart;
}
