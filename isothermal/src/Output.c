#include "jupiter.h"

void WriteCommSource (Field, filename)
     ScalarField *Field;
     char *filename;
{
  FILE *hdl;
  long gncell[3], Size, i, j, k, mg, m, stride[3], ex[3], ngh[3];
  real *buffer;
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = (WriteGhosts ? Nghost[i] : 0);
    ex[i] = Nghost[i]-ngh[i];
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
  }
  Size = (gncell[0]-2*ex[0])*(gncell[1]-2*ex[1])*(gncell[2]-2*ex[2]);
  m = 0;
  buffer = prs_malloc (Size*sizeof(real));
  for (k = ex[2]; k < gncell[2]-ex[2]; k++) {
    for (j = ex[1]; j < gncell[1]-ex[1]; j++) {
      for (i = ex[0]; i < gncell[0]-ex[0]; i++) {
	mg = i*stride[0]+j*stride[1]+k*stride[2];
	buffer[m] = (real)(Field->desc->SrcCommPerp[mg]);
	m++;
      }
    }
  }
  fwrite (buffer, sizeof(real), Size, hdl);
  free (buffer);
  fclose (hdl);
}

void WriteWorkArray (field, filename)
     real *field;
     char *filename;
{
  FILE *hdl;
  long gncell[3], Size, i;
  char fullname[MAXLINELENGTH];
  sprintf (fullname, "%s%s%ld_%ld_%ld.dat", CurrentFluidPatch->Fluid->Name,\
	   filename, 0L, CurrentFluidPatch->desc->number, CurrentFluidPatch->desc->level);
  printf ("Dumping %s..", fullname);
  fflush (stdout);
  hdl = prs_opend (fullname);
  for (i = 0; i < 3; i++)
    gncell[i] = CurrentFluidPatch->desc->gncell[i];
  Size = gncell[0] * gncell[1] * gncell[2];
  fwrite (field, sizeof(real), Size, hdl);
  fclose (hdl);
  printf ("done\n");
  fflush (stdout);
}

void WriteCommDest (Field, filename)
     ScalarField *Field;
     char *filename;
{
  FILE *hdl;
  long gncell[3], Size, i, j, k, mg, m, stride[3], ex[3], ngh[3];
  real *buffer;
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = (WriteGhosts ? Nghost[i] : 0);
    ex[i] = Nghost[i]-ngh[i];
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
  }
  Size = (gncell[0]-2*ex[0])*(gncell[1]-2*ex[1])*(gncell[2]-2*ex[2]);
  m = 0;
  buffer = prs_malloc (Size*sizeof(real));
  for (k = ex[2]; k < gncell[2]-ex[2]; k++) {
    for (j = ex[1]; j < gncell[1]-ex[1]; j++) {
      for (i = ex[0]; i < gncell[0]-ex[0]; i++) {
	mg = i*stride[0]+j*stride[1]+k*stride[2];
	buffer[m] = (real)(Field->desc->DestCommPerp[mg]);
	m++;
      }
    }
  }
  fwrite (buffer, sizeof(real), Size, hdl);
  free (buffer);
  fclose (hdl);
}

void WriteScalarField (Field, filename)
     ScalarField *Field;
     char *filename;
{
  FILE *hdl;
  long gncell[3], Size, i, j, k, mg, m, stride[3], ex[3], ngh[3];
  real *buffer;
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = (WriteGhosts ? Nghost[i] : 0);
    ex[i] = Nghost[i]-ngh[i];
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
  }
  Size = (gncell[0]-2*ex[0])*(gncell[1]-2*ex[1])*(gncell[2]-2*ex[2]);
  m = 0;
  buffer = prs_malloc (Size*sizeof(real));
  for (k = ex[2]; k < gncell[2]-ex[2]; k++) {
    for (j = ex[1]; j < gncell[1]-ex[1]; j++) {
      for (i = ex[0]; i < gncell[0]-ex[0]; i++) {
	mg = i*stride[0]+j*stride[1]+k*stride[2];
	buffer[m] = Field->Field[mg];
	m++;
      }
    }
  }
  fwrite (buffer, sizeof(real), Size, hdl);
  free (buffer);
  fclose (hdl);
}

void WriteVectorField (Field, filename)
     VectorField *Field;
     char *filename;
{
  FILE *hdl;
  long gncell[3], stride[3], Size, d, i, j, k, m, mg, ex[3], ngh[3];
  real *buffer;
  hdl = prs_opend (filename);
  for (i = 0; i < 3; i++) {
    ngh[i] = (WriteGhosts ? Nghost[i] : 0);
    ex[i] = Nghost[i]-ngh[i];
  }
  for (i = 0; i < 3; i++) {
    gncell[i] = Field->desc->gncell[i];
    stride[i] = Field->desc->stride[i];
  }
  Size = (gncell[0]-2*ex[0])*(gncell[1]-2*ex[1])*(gncell[2]-2*ex[2]);
  buffer = prs_malloc (Size*sizeof(real));
/*printf("filename %s \n",filename);
printf("strides: %ld %ld %ld \n",stride[0],stride[1],stride[2]);
printf("ex,gncell-ex: %ld %ld  \n",ex[0],gncell[0]-ex[0]);
printf("ex,gncell-ex: %ld %ld  \n",ex[1],gncell[1]-ex[1]);
printf("ex,gncell-ex: %ld %ld  \n",ex[2],gncell[2]-ex[2]);
*/
  for (d = 0; d < NDIM; d++) {
    m = 0;
    for (k = ex[2]; k < gncell[2]-ex[2]; k++) {
      for (j = ex[1]; j < gncell[1]-ex[1]; j++) {
	for (i = ex[0]; i < gncell[0]-ex[0]; i++) {
	  mg = i*stride[0]+j*stride[1]+k*stride[2];
	  buffer[m] = Field->Field[d][mg];
	  m++;
	}
      }
    }
    fwrite (buffer, sizeof(real), Size, hdl);
  }
  free (buffer);
  fclose (hdl);
}

void WriteFluid (patch, number)
     FluidPatch *patch;
     long number;
{
  char filename[MAXLINELENGTH];
  long i;
  ScalarField *sf=NULL;
  VectorField *vf=NULL;
  boolean WriteScalar=YES;
  fflush (stdout);
  /*   The following lines should be uncommented for debugging purposes only
    sprintf (filename, "gascomm%ld_%ld_%ld.dat", number,\
    patch->desc->number, patch->desc->level);
    WriteCommSource (patch->Density, filename); 
    sprintf (filename, "gasdest%ld_%ld_%ld.dat", number,\
    patch->desc->number, patch->desc->level);
    WriteCommDest (patch->Density, filename); 
  */
  for (i = 0; i < 4; i++) {
    switch (i) {
    case 0: 
      sf = patch->Density;
      WriteScalar = TRUE;
      break;
    case 1: 
      sf = patch->Energy;
      WriteScalar = TRUE;
      break;
    case 3: 			/* debugging purpose only */
      sf = patch->Potential;
      WriteScalar = TRUE;
      break;
    case 2: 
      vf = patch->Velocity;
      WriteScalar = FALSE;
      break;
    }
    if (WriteScalar) {
      sprintf (filename, "%s%s%ld_%ld_%ld.dat",\
	       patch->Name, sf->Name, number,\
	       sf->desc->number, sf->desc->level);
      WriteScalarField (sf, filename); 
    }
    else {
      sprintf (filename, "%s%s%ld_%ld_%ld.dat",\
	       patch->Name, vf->Name, number,\
	       vf->desc->number, vf->desc->level);
      WriteVectorField (vf, filename); 
    }
  }
}

void MakeDir (number)
     long number;
{
  char command[MAXLINELENGTH];
  setout (number);
  if (!CPU_Rank || AllCPUs) {
    sprintf (command, "mkdir -p %s", OutputDir);
    system (command);
  }
  MPI_Barrier (MPI_COMM_WORLD);
}

void Write (number)
     long number;
{
  tGrid_CPU *grid;
  FluidPatch *Fluid;
  prs_msg ("Output (dir: %s)  #%ld/%ld...", OUTPUTDIR, number, NTOT/NINTERM);
  pInfo ("Output (dir: %s) #%ld/%ld...", OUTPUTDIR, number, NTOT/NINTERM);
  WriteDescriptor (number);
  grid = Grid_CPU_list;
  while (grid != NULL) {
    if (grid->cpu == CPU_Rank) {
      Fluid = grid->Fluid;
      while (Fluid != NULL) {
	WriteFluid (Fluid, number);
	Fluid = Fluid->next;
      }
    }
    grid=grid->next;
  }
  prs_msg ("Done\n");
  pInfo ("Done\n");
}
