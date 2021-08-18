#include "jupiter.h"

void MonitorConservative () {
  real Mass, Mom[3]={0.,0.,0.};
  char command[MAXLINELENGTH];
  long i;
  FILE *log;
  Mass = TotalMass ();
  for (i = 0; i < NDIM; i++)
    Mom[i] = TotalMomentum (i);
  if (!CPU_Rank) {
    if ( CurrentOutputNumber == 0 ) { // a new file is created for each simulation
      sprintf(command, "cd %s; mv -f conslog.dat conslog.dat.old", OUTPUTDIR);
      system (command);
      log = prs_open ("conslog.dat");
      }
    else {
      log = prs_opena ("conslog.dat");
      }
    fprintf (log, "%ld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",\
	     CurrentOutputNumber, GlobalDate, Mass,\
	     Mom[0], Mom[1], Mom[2]);
    fclose (log);
  }
  pInfo ("Mass     total   = %.18g\n", Mass);
  for (i = 0; i < NDIM; i++)
    pInfo ("Momentum total %ld = %.18g\n", i, Mom[i]);
}

void MonitorKineticEnergy () {
  real ekin;
  char command[MAXLINELENGTH];
  FILE *log;
  ekin = TotalKineticEnergy ();
  if (!CPU_Rank) {
    if ( CurrentOutputNumber == 0 ) { // a new file is created for each simulation
      sprintf(command, "cd %s; mv -f ekinlog.dat ekinlog.dat.old", OUTPUTDIR);
      system (command);
      log = prs_open ("ekinlog.dat");
    }
    else {
      log = prs_opena ("ekinlog.dat");
    }
    fprintf (log, "%ld\t%.17g\t%.17g\n",\
	     CurrentOutputNumber, GlobalDate, ekin);
    fclose (log);
  }
  pInfo ("Total kinetic energy = %.18g\n", ekin);
}

/*
 monitoring the torque exerted by the planet on the disk (in torque.dat)
 */
void MonitorTorque () {
  int Nsides = 2;
  int Nlimits = 3;
  int sides[] = {-1,+1};				// sides[i] = inner torque, outer torque, global torque
  real limits[] = {0., 1.0, 0.5};	// planet_distance >= limits[j] * Hills_radius
  real *totaltorques;
  FILE *log;
  char command[MAXLINELENGTH];
  int i,j;
  totaltorques = (real *)prs_malloc(Nsides*Nlimits*sizeof(real));
  if ( CurrentOutputNumber != 0 ) { // the potential is undefined at output 0
    TotalTorque(totaltorques,Nsides,sides,Nlimits,limits,1.); // 1. is a dummy value
  } else {
    for (j=0; j<Nlimits; j++) {
      for (i=0; i<Nsides; i++) {
	totaltorques[i*Nlimits+j] = 0.; // the first line of torque.dat is fake
      }
    }
  }
  if (!CPU_Rank) {
    if ( CurrentOutputNumber == 0 ) { // a new file is created for each simulation
      sprintf(command, "cd %s; mv -f torque.dat torque.dat.old", OUTPUTDIR);
      system (command);
      log = prs_open ("torque.dat");
    }
    else {
      log = prs_opena ("torque.dat");
    }
    fprintf (log, "%.17g", GlobalDate);
    for (j=0; j<Nlimits; j++) {
      for (i=0; i<Nsides; i++) {
	fprintf(log,"\t%.17g",totaltorques[i*Nlimits+j]);
      }
    }
    fprintf(log,"\n");
    fclose (log);
  }
  free (totaltorques);
}

void MonitorTorqueZ () {
  int Nz;
  real *coords;
  real minlev_size=1e20, MinCellSize=1e20, maxlev_size, maxlev_min;
  real *totaltorques, *colatitudes;
  FILE *log;
  char command[MAXLINELENGTH];
  int i,j;
  tGrid *item;
  item = GridList;
  while (item != NULL){
    if (item->level == 0) {
      maxlev_size = item->size[_COLAT_];
      maxlev_min = item->corner_min[_COLAT_];
    }
    if (item->size[_COLAT_] < minlev_size){
      minlev_size = item->size[_COLAT_];
      //MinCellSize = minlev_size / item->ncell[_COLAT_];
      MinCellSize = item->Edges[_COLAT_][1]-item->Edges[_COLAT_][0];
    }
    item=item->next;
  }
  Nz = (int)round(maxlev_size / MinCellSize);
  totaltorques = (real *)prs_malloc(Nz*sizeof(real));
  colatitudes = (real *)prs_malloc(Nz*sizeof(real));
  for (i=0; i<Nz; i++) {
      colatitudes[i] = maxlev_min + (i+0.5)*MinCellSize; // coordinates on the first line
  }
  if ( CurrentOutputNumber != 0 ) { // the potential is undefined at output 0
    TotalTorque(totaltorques,1,1,Nz,colatitudes, MinCellSize);
  } else {
    for (i=0; i<Nz; i++) {
      totaltorques[i] = colatitudes[i]; // coordinates on the first line
    }
  }
  if (!CPU_Rank) {
    if ( CurrentOutputNumber == 0 ) { // a new file is created for each simulation
      sprintf(command, "cd %s; mv -f torqueZ.dat torqueZ.dat.old", OUTPUTDIR);
      system (command);
      log = prs_open ("torqueZ.dat");
      fprintf(log,"%i\n",Nz+1);
    }
    else {
      log = prs_opena ("torqueZ.dat");
    }
    fprintf (log, "%.17g", GlobalDate);
    for (i=0; i<Nz; i++) {
      fprintf(log,"\t%.17g",totaltorques[i]);
    }
    fprintf(log,"\n");
    fclose (log);
  }
  free(totaltorques);
  free(colatitudes);
}
