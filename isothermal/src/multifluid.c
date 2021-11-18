/** \file multifluid.c

Performs preliminary work about multifluid initialization.

*/

#include "jupiter.h"

void MultiFluid ()
{
  char *Fluids, *InitCodes, *InitCodesEq;
  long i, NbInit, NbInitEq;
  char *c; 
  boolean again=NO;
  int ngas, ndust = 0;
  char test;
    
  Fluids = prs_malloc (sizeof(char)*MAXLINELENGTH);
  
  NbFluids = strlen(FLUIDS); //number of fluids
  strcpy (Fluids, FLUIDS); //FLUIDS is the parameter from the parameter file 
  for (i = 0; i < NbFluids; i++){
    if (Fluids[i] == 48){ // ASCII code of "0" is 48
      ++ngas;
      strcpy (FluidName[i], "gas"); 
    }else{
      if (Fluids[i] == 49){ // ASCII code of "1" is 49
      ++ndust;
        strcpy (FluidName[i], "dust",ndust);
      }else{
        printf("Unknown fluid type: %c\n",FLUIDS[i]);
      }
    } 
  }

  printf("The total number of fluids is %d, %d gas fluids and %d dust fluids  \n", NbFluids, ngas, ndust);
  for (i = 0; i < NbFluids; i++) {
    printf ("fluid %c: %s\n",Fluids[i],FluidName[i]);
  }
  
  do {
    c = strchr(Fluids, '/'); //strchr finds the first occurrence of a character in a string
    printf("c variable: %s \n",c);
    
    if (c != NULL) { // more than one fluid
      ++again;
      *c = 0; // pointer to the first letter 
    } else
      again = NO;
    //strcpy (FluidName[NbFluids++], Fluids); //name of fluid 1

    if (c != NULL)
      Fluids = c+1;
  } while (again);


  


  InitCodes = prs_malloc (sizeof(char)*MAXLINELENGTH);
  strcpy (InitCodes, INITCODE);
  NbInit=0;
  do {
    c = strchr(InitCodes, '/');
    if (c != NULL) {
      again = YES;
      *c = 0;
    } else
      again = NO;
    strcpy (InitCodeNames[NbInit++], InitCodes);
    if (c != NULL)
      InitCodes = c+1;
  } while (again);
  if (NbInit > NbFluids) pWarning ("There are more initialization codes than fluids !\n");

  InitCodesEq = prs_malloc (sizeof(char)*MAXLINELENGTH);
  strcpy (InitCodesEq, INITHYDROSTAT);
  NbInitEq=0;
  do {
    c = strchr(InitCodesEq, '/');
    if (c != NULL) {
      again = YES;
      *c = 0;
    } else
      again = NO;
    strcpy (InitCodeNamesEq[NbInitEq++], InitCodesEq);
    if (c != NULL)
      InitCodesEq = c+1;
  } while (again);
  if (NbInitEq > NbFluids) pWarning ("There are more equilibrium initialization codes than fluids !\n");
  for (i = 0; i < NbFluids; i++) {
    if (i >= NbInitEq)
      strcpy (InitCodeNamesEq[i], InitCodeNamesEq[NbInitEq-1]);
    if (i >= NbInit)
      strcpy (InitCodeNames[i], InitCodeNames[NbInit-1]);
  }
  pInfo ("%ld-fluid calculation:\n", NbFluids);
  for (i = 0; i < NbFluids; i++) {
    pInfo ("%s, init code: %s, equilibrium code: %s\n",\
	   FluidName[i], InitCodeNames[i], InitCodeNamesEq[i]);
  }
}

void SetFluidProperties (fluid)
     FluidPatch *fluid;
{
  switch (InitMode) {
  case STANDARD:
    InitCode = fluid->InitCode;
    break;
  case EQUILIBRIUM:
    InitCode = fluid->InitCodeEq;
    break;
  }
}

void FluidCoupling (item, dt)	/* A simple implicit function for 2-fluid situations */
     tGrid_CPU *item;
     real dt;
{
  real ***v, **d, **cs;
  long gncell[3], stride[3], m, i, j, k, l;
  real e, d1, d2, cs2, v1, v2, idenom;
  FluidPatch *fluid, *reffluid;

  real *_radius, *_colat;
  real radius,colat, omegakep, C;

  if (NbFluids < 2) return;
  if (item->cpu != CPU_Rank) return;
  if (NbFluids > 2) {
    pWarning ("Coupling of more than two fluids not implemented.\n");
    return;
  }
  /*e = COUPLING*dt;*/
  v = (real ***)prs_malloc (sizeof(real *) * NbFluids);
  d = (real **)prs_malloc (sizeof(real **) * NbFluids);
  cs = (real **)prs_malloc (sizeof(real **) * NbFluids);
  fluid = item->Fluid;
  /*reffluid = item->Fluid;*/
  i = 0;
  while (fluid != NULL) {
    d[i] = fluid->Density->Field;
    v[i] = fluid->Velocity->Field;
    cs[i]= fluid->Energy->Field;
    fluid = fluid->next;
    i++;
  }
  getgridsize (item, gncell, stride);
  _radius = item->Fluid->desc->Center[_RAD_];
  _colat = item->Fluid->desc->Center[_COLAT_];
  for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	d1 = d[0][m];  /* dust density*/
	d2 = d[1][m]; /* gas density*/

  cs2 = cs[1][m]; /* square of the local sound speed of the gas*/
  radius = _radius[m]; /* radius coordinate */
  colat = _colat[m];
  omegakep = OMEGAFRAME*sin(colat)/(sqrt(radius)*sqrt(radius)*sqrt(radius)); /* local keplerian frequency */

  /*printf("%lg \n",radius);
  printf("%lg \n",radius2);
  printf(".... \n");/*
	/*idenom = 1./(1.+(d1+d2)*e);*/

  /*-----------------------------------------------------------------*/
  /*constant coupling constant*/
   /*C = COUPLING;*/
   /*-----------------------------------------------------------------*/

 /*-----------------------------------------------------------------*/
  /*constant dust particle zize in 3D*/
  //C=1.6*sqrt(cs2)/DUSTSIZE;
  /*where DUSTSIZE is the product of partcle diameter and dust solid density in code units*/
   /*-----------------------------------------------------------------*/


   /*constant dust particle zize in 2D*/
   //C=0.64*omegakep/DUSTSIZE;
   /*where DUSTSIZE is the product of partcle diameter and dust solid density in code units*/
    /*-----------------------------------------------------------------*/

  /*-----------------------------------------------------------------*/
  /*constant Stokes number*/
  // C = omegakep/(d2*STOKESNUMBER);
   /*-----------------------------------------------------------------*/

   if (NDIM==1){
     C = 1.0/STOKESNUMBER;
   }

   if (NDIM == 2){
     omegakep = OMEGAFRAME*(sqrt(radius)*sqrt(radius)*sqrt(radius));

     if(constSt==TRUE){

       C = omegakep/(d2*(STOKESNUMBER));//const stokes number
     }
     else{ //constant particle size

     C=0.64*omegakep/DUSTSIZE;
   }

     C=C*(1.0+pow((DUSTDENSFLOOR*10./d1),5)); //smooth coupling limiter
   }


 if (NDIM ==3){

   if(constSt==TRUE){
     C = omegakep/(d2*STOKESNUMBER);//const stokes number
   }
   else{ //constant particle size
     C=1.6*sqrt(cs2)/DUSTSIZE; //const dust particle size
   }

  if (d1>=DUSTDENSFLOOR){
  C=C*(1.0+pow((DUSTDENSFLOOR*100./d1),5)); //smooth coupling limiter
  }
}


	for (l = 0; l < NDIM; l++) {
	  v1 = v[0][l][m];
	  v2 = v[1][l][m];
    e = C*dt;  /*subsonic drag*/
    /* e = C*sqrt(1.0+0.221*((v1-v2)*(v1-v2))/(cs2))*dt;  supersonic drag */
    idenom = 1./(1.+(d1+d2)*e);
	  v[0][l][m] = (v1*(1.+d1*e)+v2*d2*e)*idenom;
	  //v[1][l][m] = (v2*(1.+d2*e)+v1*d1*e)*idenom;
    v[1][l][m] = v[1][l][m]; // no feedback onto the gas





    if (d1<=DUSTDENSFLOOR){  // hard coupling limiter
      v[0][l][m]=v2;
      v[1][l][m]=v2;
    }

    real vellim = 3.0;

    if (v[0][l][m]>vellim){
      v[0][l][m]=vellim;
    }
    if (v[0][l][m]<-vellim){
      v[0][l][m]=-vellim;
    }


    if (v[1][l][m]>vellim){
      v[1][l][m]=vellim;
    }
    if (v[1][l][m]<-vellim){
      v[1][l][m]=-vellim;
    }




	      }
      }
    }
  }
}
