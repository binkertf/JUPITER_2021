#include "jupiter.h"
#include <time.h>

void prs_exit (numb)
     int numb;
{
  MPI_Finalize ();
  exit (numb);
}

void prs_error (const char *template, ...)
{
  va_list ap;
  if (!CPU_Rank) {
    va_start (ap, template);
    vfprintf (stderr, template, ap);
    fprintf (stderr, "\nAborted.\n");
    va_end (ap);
    fflush (stderr);
  }
  prs_exit (EXIT_FAILURE);
}

void prs_end (const char *template, ...)
{
  va_list ap;
  pInfo (template);
  if (!CPU_Rank) {
    va_start (ap, template);
    vfprintf (stdout, template, ap);
    fprintf (stdout, "\n");
    va_end (ap);
    fflush (stderr);
  }
  prs_exit (EXIT_SUCCESS);
}

void prs_stderr (const char *template, ...)
{
  va_list ap;
  //  pError (template); Pb at start if no parameter file... 
  if (!CPU_Rank) {
    va_start (ap, template);
    vfprintf (stderr, template, ap);
    va_end (ap);
    fflush (stderr);
  }
}

void *prs_malloc (number_of_bytes)
     size_t number_of_bytes;
{
  void *ptr;
  long i;
  if (number_of_bytes <= 0) 
    return NULL;
  ptr = malloc (number_of_bytes);
  if (ptr == NULL) prs_error ("Not enough memory.");
  for (i = 0; i < number_of_bytes; i++)
    *((char *)ptr+i) = 0;
  return ptr;
}

real **multiple_alloc_1D (nb,size)
     long nb,size;
{
  long i;
  real **ptr;
  ptr = prs_malloc(nb*sizeof(real *));
  ptr[0] = prs_malloc(nb*size*sizeof(real));
  for (i = 1; i < nb; i++)
    ptr[i] = ptr[i-1]+size;
  return ptr;
}

FILE *prs_open (filename)
char *filename;
{
  char fullname[MAXLINELENGTH];
  char command[MAXLINELENGTH];
  char outputdir[MAXLINELENGTH];
  char *home;
  char rcfile[MAXLINELENGTH];
  FILE *handle, *rc;
  time_t tloc;
  static boolean AlreadyWrittenLastOut = NO;
  sprintf (fullname, "%s/%s", OUTPUTDIR, filename); 
/* The 'slash' is required because this function can be called while
   the code is reading the parameter file (e.g. to open the
   'warning.log' file) and at this point a trailing slash may has not
   been appended to OUTPUTDIR */
  handle = fopen (fullname, "w");
  if (handle == NULL) {
				/* we shall try again, after creating
				   the OUTPUTDIR directory */
    sprintf (command, "mkdir -p %s", OUTPUTDIR);
    system (command);
    handle = fopen (fullname, "w");
    if (handle == NULL)
      prs_error ("I am really unable to open %s.", fullname);
  }
  /* We now keep a track for convenience of the output directory */
  /* This allows to define an alias that is a shortcut to the output */
  /* directory. We prefer not to use an environment variable to do so */
  /* so that the output directory can be accessed through any session, */
  /* e.g. alias jo='cd `tail $HOME/.jupiterrc/lastout`' under bash. */
  /* We wish this path to be an absolute path, so if we detect a local */
  /* path, we concatenate it to the current working directory */
  if ((!CPU_Rank) && (!AlreadyWrittenLastOut)) {
    home = getenv ("HOME");
    strcpy(rcfile, home);
    strcat(rcfile, "/.jupiterrc");
    mkdir (rcfile, 0755);
    strcat (rcfile, "/lastout");
    rc = fopen (rcfile, "a");
    if (rc) { 			/* Silent if cannot create file */
      time (&tloc);
      fprintf (rc, "%s", ctime(&tloc));
      if (*OUTPUTDIR == '/') {
	fprintf (rc, "%s\n", OUTPUTDIR);
      } else {
	getcwd (outputdir, MAXLINELENGTH);
	strcat (outputdir, "/");
	strcat (outputdir, OUTPUTDIR);
	fprintf (rc, "%s\n", outputdir);
      }
      AlreadyWrittenLastOut = YES;
      fclose (rc);
    }
  }
  return handle;
}

FILE *prs_openi (filename)
char *filename;
{
  char fullname[MAXLINELENGTH];
   FILE *handle;
  sprintf (fullname, "%s%s", InputPath, filename); 
  handle = fopen (fullname, "w");
  if (handle == NULL) {
    prs_error ("I am unable to open %s.", fullname);
  }
  return handle;
}

FILE *prs_opena (filename)
char *filename;
{
  char fullname[MAXLINELENGTH];
  FILE *handle;
  sprintf (fullname, "%s%s", OUTPUTDIR, filename);
  handle = fopen (fullname, "a");
  if (handle == NULL)
    prs_error ("Unable to open %s.", fullname);
  return handle;
}

void setout (number)
     long number;
{
  sprintf (OutputDir, "%soutput%05ld/", OUTPUTDIR, number);
}

FILE *prs_opend (filename)
char *filename;
{
  char fullname[MAXLINELENGTH];
  FILE *handle;
  sprintf (fullname, "%s%s", OutputDir, filename);
  handle = fopen (fullname, "w");
  if (handle == NULL)
    prs_error ("Unable to open %s.", fullname);
  return handle;
}

FILE *prs_openrd (filename)
char *filename;
{
  char fullname[MAXLINELENGTH];
  FILE *handle;
  sprintf (fullname, "%s%s", OutputDir, filename);
  handle = fopen (fullname, "r");
  if (handle == NULL)
    prs_error ("Unable to open %s.", fullname);
  return handle;
}

FILE *prs_openr (filename)
char *filename;
{
  char fullname[MAXLINELENGTH];
  FILE *handle;
  sprintf (fullname, "%s%s", OUTPUTDIR, filename);
  handle = fopen (fullname, "r");
  if (handle == NULL)
    prs_error ("Unable to open %s.", fullname);
  return handle;
}

inline void getgridsize (grid, gncell, stride)
     tGrid_CPU *grid;
     long gncell[3];
     long stride[3];
{
  int k;
  for (k = 0; k < 3; k++) {	/* do not change this 3 to a NDIM. */
    gncell[k] = grid->gncell[k];
    stride[k] = grid->stride[k];
  }
}

inline void getgridsizes (grid, gncell, stride, gncells, strides)
     tGrid_CPU *grid;
     long gncell[3], gncells[3];
     long stride[3], strides[3];
{
  int k;
  for (k = 0; k < 3; k++) {	/* do not change this 3 to a NDIM. */
    gncell[k] = grid->gncell[k];
    stride[k] = grid->stride[k];
    if (HydroStaticEnforce) {
      gncells[k] = grid->Fluid->Rho_eq_c->gncell[k];
      strides[k] = grid->Fluid->Rho_eq_c->stride[k];
    }
  }
}

void memcpystride (dest, src, n, stride)
     real *src, *dest;
     long n, stride;
{
  if (stride == 1) {
    memcpy (dest, src, (size_t)(n*sizeof(real)));
  } else {
    while (n--) {
      dest[0] = src[0];
      dest++;
      src += stride;
    }
  }
}

void multarray (array, factor, n)
     real *array, factor;
     long n;
{
  while (n--) 
    *array++ *= factor;
}

void arraymult (array1, array2, n)
     real *array1, *array2;
     long n;
{
  while (n--)
    *(array1++) *= *(array2++);
}

void mpiwarning () {
#ifdef _PARALLEL
  fprintf (stdout, "Warning : parallel version.\n");
  fprintf (stdout, "It should be run only with the appropriate\n");
  fprintf (stdout, "tool, such as mpirun, prun, etc.\n");
#endif
}

void prs_msg (const char *template, ...)
{
  va_list ap;
  if (CPU_Rank) return;
  va_start (ap, template);
  vfprintf (stdout, template, ap);
  va_end (ap);
  fflush (stdout);
}

void cpTriplet (src, dest)
     real src[3], dest[3];
{
  long i;
  for (i = 0; i < 3; i++)
    dest[i] = src[i];
}

void resetTriplet (t)
     real t[3];
{
  long i;
  for (i = 0; i < 3; i++)
    t[i] = 0.0;
}

void ForAllPatches (func)
     void (*func)();
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  item = Grid_CPU_list;
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      Fluid = item->Fluid;
      while (Fluid != NULL) {
	func(Fluid);
	Fluid = Fluid->next;
      }
    }
    item = item->next;
  }
}

void SetToZero (ptr, size)
     real *ptr;
     long size;
{
  long i;
  for (i = 0; i < size; i++)
    ptr[i] = 0.0;
}

void swapl (a, b)
     long *a, *b;
{
  long c;
  c = *b;
  *b = *a;
  *a = c;
}
