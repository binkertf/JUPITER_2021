#include "jupiter.h"

void ReadField (fp, NbRestart)
     FluidPatch *fp;
     long NbRestart;
{
  char filename[MAXLINELENGTH];
  char varname[MAXLINELENGTH];
  long nvar, ndim=0, i, j, k, n, size[3], Size[3], dim, m;
  unsigned long q;
  long pcmin[3], pcmax[3], gncell[3], stride[3], gsize;
  int shakehand;
  MPI_Status stat;
  FILE *in;
  real *field=NULL, *buffer;
  size_t sz, exp_sz, size_ratio;
  nvar = 3;
  if (Stellar) nvar=9;
  getgridsize (fp->desc, gncell, stride);
  for (i = 0; i < 3; i++) {
    pcmin[i] = fp->desc->pcorner_min[i];
    pcmax[i] = fp->desc->pcorner_max[i];
    size[i] = fp->desc->ncell[i];
    Size[i] = fp->desc->Parent->ncell[i];
  }    
  gsize = gncell[0]*gncell[1]*gncell[2];
  buffer = prs_malloc(sizeof(real)*Size[0]);
  for (i = 0; i < nvar; i++) {
    switch (i) {
    case 0:
      field = fp->Density->Field;
      strcpy (varname, fp->Density->Name);
      ndim = 1;
      break;
    case 1:
      field = fp->Energy->Field;
      strcpy (varname, fp->Energy->Name);
      ndim = 1;
      break;
    case 2:
      field = fp->Velocity->Field[0];
      strcpy (varname, fp->Velocity->Name);
      ndim = NDIM;
      break;
    if (Stellar) {
      case 3:
        field = fp->EnergyRad->Field;
        strcpy (varname, fp->EnergyRad->Name);
        ndim = 1;
        break;
      case 4:
        field = fp->OpticalDepth->Field;
        strcpy (varname, fp->OpticalDepth->Name);
        ndim = 1;
        break;
      case 5:
        field = fp->StellarHeating->Field;
        strcpy (varname, fp->StellarHeating->Name);
        ndim = 1;
        break;
      case 6:
        field = fp->Temperature->Field;
        strcpy (varname, fp->Temperature->Name);
        ndim = 1;
        break;
      case 7:
        field = fp->Opacity->Field;
        strcpy (varname, fp->Opacity->Name);
        ndim = 1;
        break;
      case 8:
        field = fp->Gamma->Field;
        strcpy (varname, fp->Gamma->Name);
        ndim = 1;
        break;
      }
    }
    sprintf (filename, "%soutput%05ld/%s%s%ld_%ld_%ld.dat",\
	     OUTPUTDIR, NbRestart, fp->Name, varname, NbRestart,\
	     fp->desc->Parent->monoCPU, fp->desc->level);
    prs_msg ("Reading %s\n", filename);
    pInfo ("Reading %s\n", filename);
    if (CPU_Rank > 0) MPI_Recv (&shakehand, 1, MPI_INT, CPU_Rank-1, 43, MPI_COMM_WORLD, &stat);
    in = fopen (filename, "r");
    fseek (in, 0L, SEEK_END);
    sz = ftell(in);
    fseek (in, 0L, SEEK_SET);
    exp_sz = Size[0]*Size[1]*Size[2]*sizeof(real)*ndim;
    if ((exp_sz % sz == 0)  && (exp_sz != sz)) {
      printf ("Size / expected : %ld / %ld\n", (long)sz, (long)exp_sz);
      printf ("Stretching along dim 0 by factor %ld\n", (long)exp_sz/(long)sz);
    }
    size_ratio = exp_sz / sz;
    if (in == NULL)
      prs_error ("I cannot read %s. Restart impossible.", filename);
    for (dim = 0; dim < ndim; dim++) {
      for (n = 0; n < Size[1]*Size[2]; n++) {
	fread (buffer, sizeof(real), Size[0]/size_ratio, in);
	for (q = 1; q < size_ratio; q++)
	  memcpy (buffer+q*Size[0]/size_ratio, buffer, Size[0]/size_ratio*sizeof(real));
	j = n % Size[1];
	k = n / Size[1];
	if ((j >= pcmin[1]) && (k >= pcmin[2]) && (j < pcmax[1]) && (k < pcmax[2])) {
	  m = (j-pcmin[1]+Nghost[1])*stride[1]+(k-pcmin[2]+Nghost[2])*stride[2];
	  memcpy (field+m+Nghost[0]+dim*gsize, buffer+pcmin[0], size[0]*sizeof(real));
	}
      }
    }
    fclose (in);
    if (CPU_Rank < CPU_Number-1) MPI_Send (&shakehand, 1, MPI_INT, CPU_Rank+1, 43, MPI_COMM_WORLD);
  }
  free (buffer);
}
  
