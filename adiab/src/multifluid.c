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
  real e, d1, d2, cs2, acs, v1, v2, idenom, dustsz, dustsolidrho, stokes, mfp;
  FluidPatch *fluid;

  real *_radius, *_colat, *_azimuth;
  real radius, colat, azimuth, omegakep, C;

  if (NbFluids < 2) return;
  if (item->cpu != CPU_Rank) return;
  if (NbFluids > 2) {
    pWarning ("Coupling of more than two fluids not implemented.\n");
    return;
  }
  //e = COUPLING*dt;

  dustsz = DUSTSIZE / R0; // diameter of dust grains in code units
  dustsolidrho = DUSTSOLIDRHO / RHO0; // solid density of dust grains in code units, typically 3 g/cm^3 in physical units

  v = (real ***)prs_malloc (sizeof(real *) * NbFluids);
  d = (real **)prs_malloc (sizeof(real **) * NbFluids);
  cs = (real **)prs_malloc (sizeof(real **) * NbFluids);
  fluid = item->Fluid;
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
  _azimuth = item->Fluid->desc->Center[_AZIM_];
  for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	d1 = d[0][m]; /* dust density*/
	d2 = d[1][m]; /* gas density*/


  cs2 = cs[1][m];
  radius = _radius[m];
  azimuth = _azimuth[m];
  colat = _colat[m];





    /*-----------------------------------------------------------------*/
  /*constant coupling constant*/
   //C = COUPLING;
   /*-----------------------------------------------------------------*/

   /*-----------------------------------------------------------------*/
  /*constant dust particle zize*/
  acs = sqrt(cs2/d2*GAMMA*(GAMMA-1.0)); //adiabatic sound speed
  C=1.334*acs/(dustsz*dustsolidrho);
  /*where DUSTSIZE is the product of partcle diameter and dust solid density in code units*/
   /*-----------------------------------------------------------------*/

   /*
   //compute the gas mean free path
   mfp = (3.9e-24/2.0e-15)*(XM0/(R0*R0))*1.0/d2;

   if ((9./4.*mfp)>dustsz){
     //compute the stokes number
     stokes = 1.334*(dustsz*dustsolidrho)/(d2*acs);
     printf ("Dustsize: %.12g\n",dustsz*R0);
     printf ("gasdensity: %.12g\n",d2*RHO0);
     printf ("mfp: %.12g\n",mfp*R0);
     printf ("Stokes: %.12g\n",stokes);
     printf("##########");

   }*/




    /*-----------------------------------------------------------------*/
    /*constant Stokes number*/
    //omegakep = OMEGAFRAME/(sqrt(radius)*sqrt(radius)*sqrt(radius)); /* local keplerian frequency */
   //C = 2.507*omegakep/(d2*STOKESNUMBER);
   /*-----------------------------------------------------------------*/


    //C=C*(1.0+1.0e-20/d1);
    //C=C*(1.0+(1.0e-23/d1)*(1.0e-23/d1)*(1.0e-23/d1)*(1.0e-23/d1)*(1.0e-23/d1)); //smooth coupling limiter
    C=C*(1.0+pow((DUSTDENSFLOOR*100/d1),5)); //smooth coupling limiter


	for (l = 0; l < NDIM; l++) {
	  v1 = v[0][l][m];
	  v2 = v[1][l][m];
    e = C*dt;
    idenom = 1./(1.+(d1+d2)*e);
	  v[0][l][m] = (v1*(1.+d1*e)+v2*d2*e)*idenom; //dust velocity
	  v[1][l][m] = (v2*(1.+d2*e)+v1*d1*e)*idenom; //gas velocity



    if (d1<=DUSTDENSFLOOR){  // hard coupling limiter
      v[0][l][m]=v2;
      v[1][l][m]=v2;
    }


    if((colat/M_PI*180.0) <= (90.-STELLDEG) || (colat/M_PI*180.0) >= (90.+STELLDEG)) {
	    //only stellar irradiation for disc above the shadow
      d[0][m]=DUSTDENSFLOOR;
      v[0][l][m]=v2;
      v[1][l][m]=v2;
  	}




	      }
      }
    }
  }
}
