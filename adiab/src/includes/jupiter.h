#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/times.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <ctype.h>
#include <signal.h>
#include <stdarg.h>
#include "types.h"
#include "gridtypes.h"
#include "prototypes.h"
#include "def.h"
#ifdef __VAR_DEF
#include "var.h"
#else
#include "extvar.h"
#endif
#ifdef __SWITCH_
#include "switches.h"
#else
#include "global_switches.h"
#endif
#ifdef _PARALLEL
#include <mpi.h>
#else
#include "mpi_dummy.h"
#endif
#ifdef _TRAP_FPE
#include <signal.h>
#include <fenv.h>
#endif
