#include "jupiter.h"

#define NCPUMAX 4096L

static long ParentPos[6*NCPUMAX];
static long Number[NCPUMAX];

void merge (number)
     long number;
{
  FILE *in, *out;
  char filename_in[MAXLINELENGTH], command[MAXLINELENGTH];
  char field_name[MAXLINELENGTH];
  char filename_out[MAXLINELENGTH], line[MAXLINELENGTH];
  long foo, i, j, k, ngrid, nvar, ncpu, nb, level, ngh, size[3], bc[6], bcc[6];
  long cpugridnb;
  real rfoo, levdate;
  setout (number);
  sprintf (command, "mv %soutput%05ld/Descriptor%ld.dat %soutput%05ld/Descriptor%ld.dat.old",\
	   OUTPUTDIR, number, number, OUTPUTDIR, number, number);
  if (!CPU_Rank) system (command);
/* WARNING: only master should fprintf to Descriptor file below. */
  MPI_Barrier (MPI_COMM_WORLD);
  prs_msg ("Merging files from output %ld\n", number);
  if (!CPU_Rank) {
    sprintf (filename_in, "%soutput%05ld/Descriptor%ld.dat.old", OUTPUTDIR, number, number);
    sprintf (filename_out, "%soutput%05ld/Descriptor%ld.dat", OUTPUTDIR, number, number);
    in = fopen (filename_in, "r");
    if (in == NULL)
      prs_error ("I don't find %s. I can't restart.", filename_in);
    out= fopen (filename_out, "w");
    if (out == NULL)
      prs_error ("I can't write %s. I can't restart.", filename_out);
    prs_msg ("Processing %s\n", filename_out);
    fscanf (in,  "%lf %ld %lf %lf", &GlobalDate, &foo, &_Small_Rho, &levdate);
    fprintf(out, "%.18g %ld %.18g %.18g\n", GlobalDate, foo, _Small_Rho, levdate);
    fscanf (in,  "%ld %ld", &i, &j);
    fprintf(out, "%ld %ld\n", i, j);
    fscanf (in,  "%ld %ld %ld %ld", &ngrid, &ncpu, &ngh, &nvar);
    fprintf(out, "%ld %d %ld %ld\n", ngrid, 1, ngh, nvar);
    for (i = 0; i < nvar; i++) {
      fscanf (in, "%s", line);
      fprintf(out,"%s\n", line);
    }
    for (i = 0; i < ngrid; i++) {
      fscanf (in, "%ld %ld %ld", &nb, &level, &ncpu);
      cpugridnb = ngrid-1-nb;
      fprintf (out, "%ld %ld %d\n", nb, level, 1);
      for (j = 0; j < 3; j++)
	fscanf (in, "%ld", size+j);
      fprintf (out, "%ld %ld %ld\n", size[0], size[1], size[2]);
      fscanf (in, "%ld %ld %ld %ld %ld %ld", &bc[0],&bc[3],&bc[1],&bc[4],&bc[2],&bc[5]);
      fprintf (out, "%ld %ld %ld %ld %ld %ld \n", bc[0], bc[3], bc[1], bc[4], bc[2], bc[5]); /* The trailing space is for matching "diff" purpose */
      for (j = 0; j < 3; j++) {
	for (k = 0; k < size[j]+2*ngh+1; k++) {
	  fscanf(in, "%lg", &rfoo);
	  fprintf(out, "%.18g ", rfoo);
	}
	fprintf(out, "\n");
      }
      for (j = 0; j < 2; j++) {
	for (k = 0; k < 6; k++) {
	  fscanf (in, "%ld", &foo);
	  fprintf (out, "%ld%s", foo, (k != 5 ? " ":"\n"));
	}
      }
      for (j = 0; j < ncpu; j++) {
	for (k = 0; k < 4; k++)
	  fscanf(in, "%ld", &foo);
	fscanf(in, "%ld", &foo);
	Number[j]=foo;
	for (k = 0; k < 3; k++) {
	  fscanf(in, "%ld", &foo);
	  ParentPos[k+6*j] = foo;
	  if (foo == 0)
	    bc[k] = 999;	/* Means CPU-grid has a lower face on the edge of the grid */
	}
	for (k = 0; k < 3; k++) {
	  fscanf(in, "%ld", &foo);
	  ParentPos[k+3+6*j] = foo;
	  if (foo == size[k])
	    bc[3+k] = 999;	/* Means CPU-grid has an upper face on the edge of the grid */
	}
	fscanf (in, "%ld %ld %ld %ld %ld %ld", &bcc[0],&bcc[3],&bcc[1],&bcc[4],&bcc[2],&bcc[5]);
	for (k = 0; k < 6; k++)
	  if (bc[k] == 999) bc[k] = bcc[k];
      }
      for (k = 0; k < NbFluids; k++) {
	sprintf (field_name, "%s%s", FluidName[k], "density");
      	merge_field (field_name,  number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
/*      if (GAMMA > 1.0)*/
	sprintf (field_name, "%s%s", FluidName[k], "energy");
      	merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
	sprintf (field_name, "%s%s", FluidName[k], "potential");
      	merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
	sprintf (field_name, "%s%s", FluidName[k], "velocity");
      	merge_field (field_name, number, ncpu, Number, ParentPos, size, level, cpugridnb, NDIM);
      	if (Stellar) {
		sprintf (field_name, "%s%s", FluidName[k], "tau");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
		sprintf (field_name, "%s%s", FluidName[k], "erad");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
		sprintf (field_name, "%s%s", FluidName[k], "stheat");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
		sprintf (field_name, "%s%s", FluidName[k], "opacity");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
		sprintf (field_name, "%s%s", FluidName[k], "temperature");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
		sprintf (field_name, "%s%s", FluidName[k], "taucell");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
    sprintf (field_name, "%s%s", FluidName[k], "gamma");
		merge_field (field_name,   number, ncpu, Number, ParentPos, size, level, cpugridnb, 1L);
      	}
        if (k==0){ //only write out these lines once to the descriptor file
          fprintf (out, "%ld %ld %ld %d %ld\n", size[0], size[1], size[2], 0, cpugridnb);
        	fprintf (out, "%d %d %d %ld %ld %ld\n", 0, 0, 0, size[0], size[1], size[2]);
        	fprintf (out, "%ld %ld %ld %ld %ld %ld \n", bc[0], bc[3], bc[1], bc[4], bc[2], bc[5]);
        }

      }
    }
    fclose (in);
    fclose (out);
  }
  prs_end ("Done.");
}
