#include "jupiter.h"

void MonitorConservative () {
  real Mass, Mom[3]={0.,0.,0.}, Etot;
  char command[MAXLINELENGTH];
  long i;
  FILE *log;
  Mass = TotalMass ();
  Etot = TotalInternalEnergy() + TotalKineticEnergy();
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
    fprintf (log, "%ld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",\
	     CurrentOutputNumber, GlobalDate, Mass,\
	     Mom[0], Mom[1], Mom[2],Etot);
    fclose (log);
  }
  pInfo ("Mass     total   = %.18g\n", Mass);
  pInfo ("Energy     total   = %.18g\n", Etot);
  for (i = 0; i < NDIM; i++)
    pInfo ("Momentum total %ld = %.18g\n", i, Mom[i]);
}

void MonitorEnergy () {
  real eint, ekin, etot;
  char command[MAXLINELENGTH];
  FILE *log;
  eint = TotalInternalEnergy ();
  ekin = TotalKineticEnergy ();
  etot = eint + ekin;
  if (!CPU_Rank) {
    if ( CurrentOutputNumber == 0 ) { // a new file is created for each simulation
      sprintf(command, "cd %s; mv -f enelog.dat enelog.dat.old", OUTPUTDIR);
      system (command);
      log = prs_open ("enelog.dat");
    }
    else {
      log = prs_opena ("enelog.dat");
    }
    fprintf (log, "%ld\t%.17g\t%.17g\t%.17g\t%.17g\n",\
	     CurrentOutputNumber, GlobalDate, eint, ekin, etot);
    fclose (log);
  }
  pInfo ("Total internal energy = %.18g\n", eint);
  pInfo ("Total kinetic energy = %.18g\n", ekin);
  pInfo ("Total energy = %.18g\n", etot);
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

void MonitorTorqueZ () {// Output torqueZ.dat, containing torque density as a function of colatitude
  int Nz;
  real minlev_size=1e20, MinCellSize=1e20, maxlev_size, maxlev_min;
  real *totaltorques, *colatitudes;
  FILE *log;
  char command[MAXLINELENGTH];
  int i;
  tGrid *item;
  item = GridList;
  while (item != NULL){ // getting minimal grid & cell size
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
    colatitudes[i] = maxlev_min + (i+0.5)*MinCellSize;  // cell-centered
  }
  if ( CurrentOutputNumber != 0 ) { // the potential is undefined at output 0
    TotalTorque(totaltorques,1, 1 ,Nz, colatitudes, MinCellSize); // 1 are dummy values
  } else {
    for (i=0; i<Nz; i++) {
      totaltorques[i] = colatitudes[i]; // write coordinates on the first line of torqueZ.dat
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


void ReadRemDustMass(NbRestart) 
     long NbRestart;
{
  FILE *fp,*fp_re,*fp_renew;
  FILE *fp_fresh;
  char command[MAXLINELENGTH];
  char command2[MAXLINELENGTH];
  char fname[MAXLINELENGTH];
  long int outNB,dum;
  real date,readmass,remmass=0.;
  if (!CPU_Rank){
  sprintf(fname,"%sdustmassRemv.dat",OUTPUTDIR);
    if (FileExists(fname)){
      if ( CurrentOutputNumber == 0 ) {// a new file is created for each simulation
        sprintf(command, "cd %s; mv -f dustmassRemv.dat dustmassRemv.dat.old", OUTPUTDIR);
        system (command);
        fp = prs_open("dustmassRemv.dat");
        fclose(fp);
      } else{
        sprintf(command2, "cd %s; mv -f dustmassRemv.dat dustmassRemv.dat.restart", OUTPUTDIR);
        system (command2);
        fp_re = prs_openr("dustmassRemv.dat.restart");
        fp_renew = prs_open("dustmassRemv.dat");
        while (fscanf(fp_re,"%ld\t%lf\t%lf",&outNB,&date,&readmass)!=EOF){
          if (outNB==NbRestart){
              remmass = readmass;
              break;
          } else{
            fprintf (fp_renew, "%ld\t%.17g\t%.17g\n",outNB, date, readmass);
          }
        }
      fclose(fp_re);
      fclose(fp_renew);
      }
    }else{
     fp_fresh = prs_open("dustmassRemv.dat");
     fclose(fp_fresh);
    }
  }
if (remmass){
  MPI_Bcast(&remmass,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  DustAvgMass = remmass;
 }
}

void ReadSublDustMass(NbRestart) 
     long NbRestart;
{
  FILE *fp,*fp_re,*fp_renew;
  FILE *fp_fresh;
  char command[MAXLINELENGTH];
  char command2[MAXLINELENGTH];
  char fname[MAXLINELENGTH];
  long int outNB,dum;
  real date,readmass,sublmass=0.;
  if (!CPU_Rank){
  sprintf(fname,"%sdustSublmass.dat",OUTPUTDIR);
    if (FileExists(fname)){
      if ( CurrentOutputNumber == 0 ) {// a new file is created for each simulation
        sprintf(command, "cd %s; mv -f dustSublmass.dat dustSublmass.dat.old", OUTPUTDIR);
        system (command);
        fp = prs_open("dustSublmass.dat");
        fclose(fp);
      } else{
        sprintf(command2, "cd %s; mv -f dustSublmass.dat dustSublmass.dat.restart", OUTPUTDIR);
        system (command2);
        fp_re = prs_openr("dustSublmass.dat.restart");
        fp_renew = prs_open("dustSublmass.dat");
        while (fscanf(fp_re,"%ld\t%lf\t%lf",&outNB,&date,&readmass)!=EOF){
          if (outNB==NbRestart){
              sublmass = readmass;
              break;
          } else{
            fprintf (fp_renew, "%ld\t%.17g\t%.17g\n",outNB, date, readmass);
          }
        }
      fclose(fp_re);
      fclose(fp_renew);
      }
    }else{
     fp_fresh = prs_open("dustSublmass.dat");
     fclose(fp_fresh);
    }
  }
if (sublmass){
  MPI_Bcast(&sublmass,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  DustSublMass = sublmass;
 }
}

void MonitorRemovedDustMass(){
  FILE *log;
 if (!CPU_Rank) {    
    log = prs_opena ("dustmassRemv.dat");
    fprintf (log, "%ld\t%.17g\t%.17g\n",CurrentOutputNumber, GlobalDate, DustAvgMass);
    fclose (log);
  }
}

void MonitorSublimatedDustMass(){
  FILE *log;
 if (!CPU_Rank) {
    log = prs_opena ("dustSublmass.dat");
    fprintf (log, "%ld\t%.17g\t%.17g\n",CurrentOutputNumber, GlobalDate, DustSublMass);
    fclose (log);
  }
}
