#include "jupiter.h"
#include <time.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/dir.h>
#include <dirent.h>
#include <sys/param.h>
#include <sys/stat.h>

/* This function is called at the beginning of a run. Its main purpose
is to check whether the user has not mistakably started with MPI a
sequential build of the code.  In that case, each processes believes
it is unique and sees only one process in the MPI_COMM_WORLD. CPU-time
wise, this is a disaster, since the CPU time consumption will be the
consumption of a sequential run multiplied by the number of processes,
not to speak about possible file access conflicts. In order to prevent
the user from running into this kind of situation, the following
detection of concommitent processes (for a sequential build) is
proposed: the process write a lock as an empty file with its PID as
its name in OUTPUTDIR/.pid.  Then, it lists this directory and
searches for files which have been written +/- 1 second at most after
or before it. If it finds such files, with names different from its
PID, it issues a warning.  In order to tidy up the .pid subdirectory,
the latter is emptied every time a parallel version of the code is run
with the corresponding OUTPUTDIR */

void checkpara () {
  char path[MAXLINELENGTH];
  char command[MAXLINELENGTH];
  int count,i;
  int seca, secb;
  long cpid, testpid;
  struct dirent **files;
  struct stat buf;
#ifdef _PARALLEL
  sprintf (command, "/bin/rm -f %s/.pid/[0-9]*[0-9]", OUTPUTDIR);
  system (command);
  return;
#endif
  sprintf (command, "mkdir -p %s/.pid", OUTPUTDIR);
  system (command);
  cpid = (long)getpid();
  sprintf (command, "touch %s/.pid/%ld", OUTPUTDIR, (long)getpid());
  system (command);
  sprintf (path, "%s/.pid", OUTPUTDIR);
  count = scandir(path, &files, fselect, alphasort);
  if (count <= 0)
    prs_error ("Internal error: no pid lock file found in %s/.pid/\nCheck 'paracheck.c'\n");
  seca=(int)(time(NULL)); 	/* Number of seconds since 1/1/1970 for current process */
  for (i = 0; i < count; i++) { /* For each file found in $OUTPUTDIR/.pid */
    sprintf (command, "%s/.pid/%s", OUTPUTDIR, files[i]->d_name); /* Check file status */
    stat (command, &buf);
    secb = (int)(buf.st_mtime);	/* Number of seconds between 1/1/1970 and file modification */
    if (abs(seca-secb) <= 1) {	/* If file is recent (we allow for one second change only) */
      testpid=atol(files[i]->d_name); /* Then we test its name */
      if (testpid != cpid) {	      /* It must be the PID of current process, otherwise */
        pWarning ("Several processes with similar starting dates");	      /* there is likely a problem */
	pWarning ("were found. These process all write data to %s.", OUTPUTDIR);
	pWarning ("This strongly suggests that you started with MPI");
	pWarning ("a SEQUENTIAL (rather than MPI) version of the code.");
      }
    }
  }
}

int fselect(entry)
  struct dirent *entry;
{
  if ((entry->d_name[0] < '0') || (entry->d_name[0] > '9'))
    return (FALSE);
  else
    return (TRUE);
}
