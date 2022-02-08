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
  real e, d1, d2, cs2, acs, acs2, acs2_2, v1, v2, idenom, dustsz, dustsolidrho, stokes, mfp, diff_f, delta, tau_s, omega_coll, lamba_mfp;
  FluidPatch *fluid;

  real *_radius, *_colat, *_azimuth;
  real radius, colat, azimuth, omegakep, C;

  if (NbFluids < 2) return;
  if (item->cpu != CPU_Rank) return;
  if (NbFluids > 2) {
    pWarning ("Coupling of more than two fluids not implemented.\n");
    return;
  }

  dustsz = DUSTSIZE / R0; // radius of dust grains in code units
  dustsolidrho = DUSTSOLIDRHO / RHO0; // solid density of dust grains in code units, typically 3 g/cm^3 in physical units
  omega_coll = 2.0 * M_PI * dustsz * dustsz; //collision cross section

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

        tau_s = 0.0;// stopping time
        
        if (Isothermal){
          switch (NDIM){
            case 1: // NDIM == 1
              omegakep = 1.0;
              if(constSt==YES){
                C = 1.0/STOKESNUMBER;
                tau_s = STOKESNUMBER / omegakep;
              }else{ //constant particle size
                prs_error ("ERROR: Constant particle size not implemented in 1-D setup. \n");
              }
              break;

            case 2: // NDIM == 2
              omegakep = 1.0/sqrt(radius*radius*radius);
              if(constSt==YES){
                C = omegakep/(d2*STOKESNUMBER);//const stokes number
                tau_s = STOKESNUMBER / omegakep;


              }else{ //constant particle size
                C=0.64*omegakep/(dustsz * dustsolidrho);
                tau_s = M_PI / 2.0 * dustsz * dustsolidrho / d2 / omegakep;
              }
              if (DIFFMODE != 1)
                C=C*(1.0+pow((DUSTDENSFLOOR*10./d1),5)); //smooth coupling limiter
              break;

            case 3: // NDIM == 3
              omegakep = 1.0*sin(colat)/(sqrt(radius)*sqrt(radius)*sqrt(radius));
              if(constSt==YES){
                tau_s = STOKESNUMBER/omegakep;
                C = omegakep/(d2*STOKESNUMBER);//const Stokes number
              }else{
                C=1.6*sqrt(cs2)/(dustsz * dustsolidrho); //const dust particle size
                tau_s = sqrt(M_PI / 8.0) * dustsz * dustsolidrho / sqrt(cs2) / d2 / omegakep;
                if (DIFFMODE !=1)
                  C=C*(1.0+pow((DUSTDENSFLOOR*100./d1),5)); //smooth coupling limiter
              }
              break;

            default:
              prs_error ("ERROR: Invalid number of dimensions.");
          }
        }else{ //radiative 
          if(constSt==YES){//constant Stokes number
            prs_error ("ERROR: Constant Stokes number not implemented in radiative setup. Use constant particle size instead -> CONSTSTOKES FALSE \n");
          }else{ //constant dust particle size
            omegakep = 1.0*sin(colat)/(sqrt(radius)*sqrt(radius)*sqrt(radius));
            acs = sqrt(cs2/d2*GAMMA*(GAMMA-1.0)); //adiabatic sound speed
            C=1.334*acs/(dustsz*dustsolidrho);
            tau_s = sqrt(GAMMA * M_PI / 8.0) * dustsz * dustsolidrho / acs / d2 / omegakep;
            //C=C*(1.0+pow((DUSTDENSFLOOR*100/d1),5)); //smooth coupling limiter
          }
        }

       if (dt>tau_s)prs_error ("\n\nWARNING: Gas-dust coupling time is not resolved.\n\
This could lead to errors. If you wish to ignore, \n\
comment out this warning in multifluid.c and \n\
PROCEED AT YOUR OWN RISK! \n\n\
t_stop = %.2e\n\
dt =     %.2e\n",tau_s,dt);

        if (constSt!=YES){ 
          lamba_mfp = 1.0 / (omega_coll* d2 / XMH); //gas mean-free-path
          if (lamba_mfp< ((4.0/9.0)*dustsz)){
            prs_error ("\n\nWARNING: Epstein might not be valid.\n\
Consider implementing Stokes drag or ignore this message at your own risk.\n");

          }
        }

        // implicit velocity update 
	      for (l = 0; l < NDIM; l++) {
	        v1 = v[0][l][m];
	        v2 = v[1][l][m];
          e = C*dt;
          idenom = 1./(1.+(d1+d2)*e);
	        v[0][l][m] = (v1*(1.+d1*e)+v2*d2*e)*idenom; //dust velocity
          if (BACKREACTION == YES){
            v[1][l][m] = (v2*(1.+d2*e)+v1*d1*e)*idenom; //gas velocity
          }



          // hard coupling limiter
          //if (d1<=DUSTDENSFLOOR){  
            //v[0][l][m]=v2;
            //v[1][l][m]=v2;
          //}


          // Here, I set the dust density to the dustfloor in the irradiated region
          if (FALSE == TRUE){ //if (!Isothermal){
            if((colat/M_PI*180.0) <= (90.-STELLDEG) || (colat/M_PI*180.0) >= (90.+STELLDEG)) {
            d[0][m]=DUSTDENSFLOOR;
            v[0][l][m]=v2;
            v[1][l][m]=v2;
  	        }
          }
	      }


        if (DIFFMODE != 1) cs[0][m] = 0.0; //set dust "energy" to zero when there is no diffusion pressure


      }
    }
  }
}


void MultifluidDiffusionPressure (item, dt)	/* Turbulent diffusion pressure in dust fluids */
     tGrid_CPU *item;
     real dt;
{
  real ***v, **d, **cs;
  long gncell[3], stride[3], m, i, j, k, l;
  real e, d1, d2, cs2, acs, acs2, acs2_min, acs2_2, v1, v2, idenom, dustsz, dustsolidrho, stokes, mfp, diff_f, delta, tau_s;
  FluidPatch *fluid;

  real *_radius, *_colat, *_azimuth;
  real radius, colat, azimuth, omegakep, C;

  if ((NbFluids < 2) || (DIFFMODE != 1) || (VISCOSITY < 1e-15)) return;
  if (item->cpu != CPU_Rank) return;
  if (NbFluids > 2) {
    pWarning ("Diffusion pressure of more than two fluids not implemented.\n");
    return;
  }
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
  
  dustsz = DUSTSIZE / R0; // diameter of dust grains in code units
  dustsolidrho = DUSTSOLIDRHO / RHO0; // solid density of dust grains in code units, typically 3 g/cm^3 in physical units

  for (i = 0; i < gncell[0]; i++) {
    for (j = 0; j < gncell[1]; j++) {
      for (k = 0; k < gncell[2]; k++) {
	      m = i*stride[0]+j*stride[1]+k*stride[2];
        d1 = d[0][m]; /* dust density*/
	      d2 = d[1][m]; /* gas density*/
        cs2 = cs[1][m];/* square of the local sound speed of the gas*/

        radius = _radius[m];
        azimuth = _azimuth[m];
        colat = _colat[m];

        if (Isothermal){//isothermal
          if (NDIM ==1){
            tau_s = STOKESNUMBER;
            cs[0][m] = VISCOSITY / tau_s; //VISCOSITY / (tau_s + VISCOSITY/cs2);
          }

          if (NDIM ==2){
          omegakep = 1.0/sqrt(radius*radius*radius);
           if(constSt==TRUE){//constant Stokes number
              delta = VISCOSITY/(sqrt(cs2)*ASPECTRATIO*radius);
              diff_f = delta/(delta+STOKESNUMBER);
              cs[0][m] = diff_f * cs2; //dust turbulent diffusion pressure

              //remove the following line after testing:
              tau_s = STOKESNUMBER / omegakep;
              cs[0][m] = VISCOSITY / tau_s;
              //until here
              
            }else{ //constant particle size
              tau_s = M_PI / 2.0 * dustsz * dustsolidrho / d2 / omegakep;
              cs[0][m] = VISCOSITY / (tau_s + VISCOSITY/cs2);
            }
          }


          if (NDIM ==3){
            if(constSt==TRUE){//constant Stokes number
              delta = VISCOSITY/(sqrt(cs2)*ASPECTRATIO*radius);
              diff_f = delta/(delta+STOKESNUMBER);
              cs[0][m] = diff_f * cs2; //dust turbulent diffusion pressure
            }else{ //constant particle size
              tau_s = sqrt(M_PI / 8.0) * dustsz * dustsolidrho / (sqrt(cs2) * d2);
              cs[0][m] = VISCOSITY / (tau_s + VISCOSITY/cs2);
            }
          }
        }else{//radiative
          acs2 = cs2/d2*(GAMMA-1.0); //adiabatic sound speed squared
          if(constSt==TRUE){
            prs_error ("ERROR: Constant Stokes number not implemented in radiative setup. Use constant particle size instead. \n");
          }else{ //constant particle size
            tau_s = sqrt(GAMMA * M_PI / 8.0) * dustsz * dustsolidrho / (sqrt(acs2) * d2);
            cs[0][m] = VISCOSITY / (tau_s + VISCOSITY/(acs2));
            //cs[0][m] = acs2; 
            //printf("%.6f \n",sqrt(cs2));
          }
        }
      }
    }
  }
}




void MultifluidDustEnergyToZero (item, dt)	/* Turbulent diffusion pressure in dust fluids */
     tGrid_CPU *item;
     real dt;
{
  real ***v, **d, **cs;
  long gncell[3], stride[3], m, i, j, k, l;
  real e, d1, d2, cs2, acs, acs2, acs2_min, acs2_2, v1, v2, idenom, dustsz, dustsolidrho, stokes, mfp, diff_f, delta, tau_s;
  FluidPatch *fluid;

  real *_radius, *_colat, *_azimuth;
  real radius, colat, azimuth, omegakep, C;

  if ((NbFluids < 2) || (DIFFMODE == 1)) return;
  if (item->cpu != CPU_Rank) return;


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
  

  for (i = 0; i < gncell[0]; i++) {
    for (j = 0; j < gncell[1]; j++) {
      for (k = 0; k < gncell[2]; k++) {
	      m = i*stride[0]+j*stride[1]+k*stride[2];

        cs[0][m] = 0.0; //set dust "energy" to zero when there is no diffusion pressure

      }
    }
  }
}