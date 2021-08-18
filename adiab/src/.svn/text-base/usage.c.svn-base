#include "jupiter.h"

void PrintUsage (execname)
     char *execname;
{
  char *end;
  end = execname+strlen(execname)-1;
  while ((end > execname) && (*(end-1) != '/')) end--;
  prs_stderr ("Usage : %s [-aceghiknptuz] [-j maxlevel] [-o option] [-S number] [-s number] [-m number] [-r \"output# grid# pos\"] [-x \"output# camera_file\"] filename\n", end);
  prs_stderr ("\n-a: All CPUs create output directories. Useful if each CPU has a different disk.\n");
  prs_stderr ("-c: monitor CFL condition or coarse ray tracing\n");
  prs_stderr ("-e: monitor kinetic energy.\n");
  prs_stderr ("-g: include ghosts when writing hydro fields to the disk\n");
  prs_stderr ("-i: superimpose mode for a restart.\n");
  prs_stderr ("-j: set maximal level of refinement.\n");
  prs_stderr ("-k: computes torque with Stockholm specifications.\n");
  prs_stderr ("-m number: merge output number (to fake mono-CPU)\n");
  prs_stderr ("-n: run disabled. Only read input files, and write first output\n");
  prs_stderr ("-o: redefines a variable (e.g. -o \"outputdir=/scratch/test/\").\n");
  prs_stderr ("-p: monitor CPU or wall clock time usage.\n");
  prs_stderr ("-r \"output# grid# position\": refine a specific grid at a given position\n");
  prs_stderr ("-s number: restart run (from mono-CPU output)\n");
  prs_stderr ("-S number: stretch run (from mono-CPU output). Not a real restart ==> output 0\n");
  prs_stderr ("-t: monitor torque.\n");
  prs_stderr ("-x \"output# camera_file\": ray tracing on given output\n");
  prs_stderr ("-z: do not monitor mass and angular momentum\n");
  prs_exit (EXIT_FAILURE);
}
