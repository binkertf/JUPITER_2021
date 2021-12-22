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
  Fluids = prs_malloc (sizeof(char)*MAXLINELENGTH);
  strcpy (Fluids, FLUIDS);
  NbFluids=0;
  do {
    c = strchr(Fluids, '/');
    if (c != NULL) {
      again = YES;
      *c = 0;
    } else
      again = NO;
    strcpy (FluidName[NbFluids++], Fluids);
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
  real radius,colat, omegakep, C, delta, diff_f, tau_s;

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
  omegakep = 1.0*sin(colat)/(sqrt(radius)*sqrt(radius)*sqrt(radius)); //OMEGAFRAME*sin(colat)/(sqrt(radius)*sqrt(radius)*sqrt(radius)); /* local keplerian frequency */

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
     omegakep = 1.0/(sqrt(radius)*sqrt(radius)*sqrt(radius)); //OMEGAFRAME/(sqrt(radius)*sqrt(radius)*sqrt(radius));

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


   }else{ //constant particle size
     C=1.6*sqrt(cs2)/DUSTSIZE; //const dust particle size

   }

  if (d1>=DUSTDENSFLOOR){
  C=C*(1.0+pow((DUSTDENSFLOOR*100./d1),5)); //smooth coupling limiter
  }
}


	for (l = 0; l < NDIM; l++) {
	  v1 = v[0][l][m];
	  v2 = v[1][l][m];
    //e = C*dt;  /*subsonic drag*/
    e = C*sqrt(1.0+0.221*((v1-v2)*(v1-v2))/(cs2))*dt;  //supersonic drag
    idenom = 1./(1.+(d1+d2)*e);
	  v[0][l][m] = (v1*(1.+d1*e)+v2*d2*e)*idenom;
	  v[1][l][m] = (v2*(1.+d2*e)+v1*d1*e)*idenom;
    //v[1][l][m] = v[1][l][m]; // no feedback onto the gas





    if (d1<=DUSTDENSFLOOR){  // hard coupling limiter
      v[0][l][m]=v2;
      v[1][l][m]=v2;
    }

	      }

        for (l = 0; l < NDIM; l++) {
        real vellim = 2.0;

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

  if (DIFFMODE == 2){
  //diffpressure
  for (i = 0; i < gncell[0]; i++) {
    for (j = 0; j < gncell[1]; j++) {
      for (k = 0; k < gncell[2]; k++) {
	       m = i*stride[0]+j*stride[1]+k*stride[2];
	        d1 = d[0][m];  /* dust density*/
	        d2 = d[1][m]; /* gas density*/
          cs2 = cs[1][m]; /* square of the local sound speed of the gas*/
          radius = _radius[m]; /* radius coordinate */

          if (NDIM ==3){
            if(constSt==TRUE){
              delta = VISCOSITY/(sqrt(cs2)*ASPECTRATIO*radius);
              diff_f = delta/(delta+STOKESNUMBER);

              cs[0][m] = diff_f * cs2; //dust diff pressure



            }else{ //constant particle size

              tau_s = sqrt(M_PI / 8.0) * DUSTSIZE / (sqrt(cs2) * d2);
              cs[0][m] = VISCOSITY / (tau_s + VISCOSITY/cs2);
            }
         }

      }
    }
  }
}





}
