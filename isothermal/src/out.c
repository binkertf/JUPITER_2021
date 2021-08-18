#include "jupiter.h"
#include <time.h>

void DumpParameters () 
{
  FILE *hdl;
  if (CPU_Rank) return;
  hdl = prs_open("params.dat");
  ListVariables (hdl);
  fclose (hdl);
  hdl = prs_open("IDL_params.dat");
  ListVariablesIDL (hdl);
  fclose (hdl);
  hdl = prs_open("coord.dat");	/* This hidden file contains 0, 1 or 2 according to the coordinate system */
  fprintf (hdl, "%ld\n", CoordType); /* (see def.h) so that IDL easily knows about it. */
  fclose (hdl);
}

void TestEndian ()
{
  FILE *hdl;
  double a;
  if (CPU_Rank) return;
  hdl = prs_open("testendian.dat");
  a = -2.3;			/* IDL will try to read this raw format file and on failure will swap endians */
  fwrite (&a, sizeof(double), 1, hdl);
  fclose (hdl);
}

void DumpSources (argc, argv)
     int argc;
     char *argv[];
{
  char commandname[MAXLINELENGTH];
  char CommandLine[MAXLINELENGTH];
  char filename[MAXLINELENGTH];
  char rcfile[MAXLINELENGTH];
  FILE *rc;
  char *ptr, *home;
  int i;
  time_t tloc;
  if (CPU_Rank) return;
  sprintf (commandname, "%s", argv[0]);
  ptr = strrchr (commandname, '/');
  if (ptr != NULL)		/* Strip any leading path information */
    strcpy (filename, ptr+1);
  else
    strcpy (filename, commandname);
  sprintf (CommandLine, "cp .src.%s.tar.bz2 %ssrc.tar.bz2", filename, OUTPUTDIR);
  system (CommandLine);
  pInfo ("***\nThis simulation was run with the ");
#ifdef _PARALLEL
  pInfo ("parallel (MPI) ");
#else
  pInfo ("sequential ");
#endif
  pInfo ("version of the JUPITER code\n");
#ifdef _PARALLEL
  pInfo ("The run was performed on %d processors\n", CPU_Number);
#endif
  pInfo ("The command line was:\n");
  home = getenv ("HOME");
  strcpy(rcfile, home);
  strcat(rcfile, "/.jupiterrc");
  mkdir (rcfile, 0755);
  strcat (rcfile, "/history");
  rc = fopen (rcfile, "a");
  if (!CPU_Rank) {
    time (&tloc);
    if (rc) fprintf (rc, "\n%s", ctime(&tloc)); 
    for (i = 0; i < argc; i++) {
      if (rc) { 	 /* Silent if cannot create file */
	fprintf (rc, "%s ", argv[i]);
	pInfo ("%s ", argv[i]);
      }
    }
  }
  if (rc) fclose (rc);
  pInfo ("\n***\n");
}
