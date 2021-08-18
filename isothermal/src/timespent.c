#include "jupiter.h"

#define NBCHRONO 10

static double InitTime[NBCHRONO], EndTime[NBCHRONO];
static char time_measured[2][128]={"Wall clock time", "CPU time"};

void ResetChronometer (long chrono)
{
  if (MonitorTimeSpent == NO) return;
  if ((chrono < 0) || (chrono >= NBCHRONO)) prs_error ("Invalid chronometer number");
  InitTime[chrono] = MPI_Wtime ();
}

void ReadChronometer (long chrono, char *msg)
{
  if (MonitorTimeSpent == NO) return;
  if ((chrono < 0) || (chrono >= NBCHRONO)) prs_error ("Invalid chronometer number");
  EndTime[chrono] = MPI_Wtime ();
#ifdef _PARALLEL
#define TIMESPENT 0
#else
#define TIMESPENT 1
#endif
  pInfo ("%s for %s: %.6g sec.\n",\
	 time_measured[TIMESPENT], msg, EndTime[chrono]-InitTime[chrono]);
}
