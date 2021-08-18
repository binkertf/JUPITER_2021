#include "jupiter.h"

//  opacity routine: Doug's Opacity Function (modified)

inline real oplin(real temp,real rho) {
  real power1,power2,power3;
  real t234,t456,t678;
  real ak1,ak2,ak3;
  real bk3,bk4,bk5,bk6,bk7;
  real ts4=0.0,rho13,rho23,ts42=0.0,ts44,ts48;
  real o3,o4,o5,o4an,o3an;
  real t2,t4,t8,t10;
  real o1,o2,o1an,o2an,o6,o7,o6an,o7an,opacity=0.0;

  //   data power1,power2,power3/4.44444444e-2, 2.381e-2, 2.267e-1/
  power1=4.44444444e-2;
  power2=2.381e-2;
  power3=2.267e-1;

  //    data t234,t456,t678/1.6e3, 5.7e3, 2.28e6/
  t234=1.6e3;
  t456=5.7e3;
  t678=2.28e6;

  //  coefficients for opacity laws 1, 2, and 3 in cgs units.
  //  data ak1,ak2,ak3/2.e-4, 2.e16, 5.e-3/
  ak1=2.e-4;
  ak2=2.e16;
  ak3=5.e-3;


  /*   coefficients for opacity laws 3, 4, 5, 6, 7, and 8 in T_4 units.
       data bk3,bk4,bk5,bk6,bk7   50., 2.e-2, 2.e4, 1.d4, 1.5d10 */

  bk3 = 50.;
  bk4 = 2.e-2;
  bk5 = 2.e4;
  bk6 = 1.e4;
  bk7 = 1.5e10;


  /*  test T against (T_23 * T_34 * T_34)**0.333333333 */

  if( temp < t234*pow(rho,power1)){
    //   different powers of temperature
    t2=temp*temp;
    t4=t2*t2;
    t8=t4*t4;
    t10=t8*t2;
    //   disjoint opacity laws
    o1=ak1*t2;
    o2=ak2*temp/t8;
    o3=ak3*temp;
    //   parameters used for smoothing
    o1an=o1*o1;
    o2an=o2*o2;
    //   smoothed and continuous opacity law for regions 1, 2, and 3.
    opacity=pow(pow(o1an*o2an/(o1an+o2an),2.0)+pow(o3/(1.0+1.e22/t10),4.0),0.25);
  }
  
  if (temp > t234*pow(rho,power1) && temp < t456*pow(rho,power2)){

    //   to avoid overflow
    ts4=1.e-4*temp;
    rho13=pow(rho,0.333333333);
    rho23=rho13*rho13;
    ts42=ts4*ts4;
    ts44=ts42*ts42;
    ts48=ts44*ts44;

    //   disjoint opacity laws for 3, 4, and 5.
    o3=bk3*ts4;
    o4=bk4*rho23/(ts48*ts4);
    o5=bk5*rho23*ts42*ts4;
    //   parameters used for smoothing
    o4an=pow(o4,4.0);
    o3an=pow(o3,4.0);
    //  smoothed and continuous opacity law for regions 3, 4, and 5.
    opacity=pow((o4an*o3an/(o4an+o3an))+pow(o5/(1.0+6.561e-5/ts48),4.0),0.25);
  }
        
  if (temp > t456*pow(rho,power2)){
    if (temp < t678*pow(rho,power3) || rho <= 1.0e-10){
      //  to avoid overflow
      ts4=1.e-4*temp;
      rho13=pow(rho,0.333333333);
      rho23=rho13*rho13;
      ts42=ts4*ts4;
      ts44=ts42*ts42;
      ts48=ts44*ts44;
      
      //   disjoint opacity laws for 5, 6, and 7.
      o5=bk5*rho23*ts42*ts4;
      o6=bk6*rho13*ts44*ts44*ts42;
      o7=bk7*rho/(ts42*sqrt(ts4));
      //   parameters used for smoothing
      o6an=o6*o6;
      o7an=o7*o7;
      //   smoothed and continuous opacity law for regions 5, 6, and 7.
      opacity=pow(pow(o6an*o7an/(o6an+o7an),2.0)+pow(o5/(1.0+pow(ts4/(1.1*pow(rho,0.04762)),10.0)),4.0),0.25);
    }
    else {
      //  no scattering!
      opacity=bk7*rho/(ts42*sqrt(ts4));
    }
  }
  return opacity;
}

//     Calculates kappa in Code units  
real kappa(real temp,real rho) {
  real xxkappa,kappa=0.0;
  int OPACITYLAW = 0; //Law selector;
  // The default value OPACITYLAW=0 is for  Bell and Lin 1994
  if(OPACITYLAW == 0){
    kappa = oplin(temp*TEMP0,rho*RHO0);
   }
   else if(OPACITYLAW == 1){
     // constant opacity
     kappa = 1.0;
   }
  // Opacity limiter
 /*  if(kappa <= 0.01)
     kappa = 0.01;*/ // This was Bertram Bitsch "fix" to avoid to low opacities for the upper layer of the disk, when there is stellar irradtion. However, it just cause harm, if one does not have stellar irradiation, because you limit the minimal opacity! So use it only on your own risk.
   xxkappa=kappa*RHO0*R0+1.0e-300;
   return xxkappa;
}

void ComputeTemperatureField() {
  long i,j,h,l,nr,ns,ni;
  real *dens, *energ, *temp,*energrad,*colat;
  real lowcol=RANGE3LOW, highcol=RANGE3HIGH;
  long stride[3],gncell[3];
  real *center[3];
  real dr,daz,dco;
  FluidWork *fw;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  for (i = 0; i < 3; i++)	/* 3, not NDIM */
    center[i] = fw->desc->Center[i];
  nr = gncell[1];
  ns = gncell[0];
  ni = gncell[2];

  dens = fw->Density;
  energ = fw->Energy;  
  temp = fw->Temperature;
  energrad = fw->EnergyRad;
  colat = fw->desc->Center[_COLAT_];

// modification made from here
	    dr=center[_RAD_][stride[1]]-center[_RAD_][0];
	    daz=center[_AZIM_][1]-center[_AZIM_][0];
	    dco=center[_COLAT_][stride[2]]-center[_COLAT_][0];
// till here

  if (HalfDisk) highcol =  1e30; //Don't set the ghost value to high
				 //altitude values if you have only
				 //half a disk (up to the equator)
  
  for( h = 0; h < ni; h++) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < ns; j++ ) {
	l = j+i*stride[1]+h*stride[2];
	// We want to set low temperature only in the ghost cells in
	// colatitude (in order to coold the disk), so we check for
	// that:
	if(colat[l] < lowcol || colat[l] > highcol ){
	  temp[l] = BOUNDTEMP/TEMP0;
	  energ[l] = CV*dens[l]*temp[l]; // Is this needed ?
	  if(Stellar) 
	    energrad[l]= BoundEcode;
	} else {
	  temp[l] = energ[l]/(CV*dens[l]);
	}
// Modification made from here -- setting up low temperatures on the planet
    	    if(center[_AZIM_][l] < 0.0+daz*2.2 && center[_AZIM_][l] > 0.0-(daz*2.2) && center[_RAD_][l] < 1.0+dr*2.2 &&  center[_RAD_][l] > 1.0-(dr*2.2) && center[_COLAT_][l] > (M_PI/2.)-(dco*2.2) &&  center[_COLAT_][l] < (M_PI/2.)) {
	      if (temp[l]*TEMP0 > 8000.0 ){
 	      temp[l]=8000.0/TEMP0;
	      energ[l] = CV*dens[l]*temp[l];
	      }
	    }
// till here
      }
    }
  }
}


void ComputeOpacity() {
  FluidWork *fw;
  long size[3], gncell[3], stride[3];
  long i, j,h, l,nr, ns, ni;
  real *opar, *opap, *opas, *temper, *dens;
  real opacity, opamin;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  
  for (i = 0; i < 3; i++)
    size[i] = fw->desc->gncell[i];
  
  nr = size[1]; // if coordpermut=213 (first coord is azimuth, 2nd rad, 3rd colat)
  ns = size[0]; // azimuth
  ni = size[2]; // co-lat 
  
  opar = fw->OpaR;
  opap = fw->OpaP;
  opas = fw->OpaS;
  temper= fw->Temperature;
  dens = fw->Density;

  for(h =0; h< ni; h++) {
    for (i = 0; i < nr;  i++) {
      for (j = 0; j < ns; j++) {
	l = j+i*stride[1]+h*stride[2];
	//	multiplication by a factor 100 of DUSTTOGAS as opacity is
	//      already calculated for a dust to gas ratio of 0.01
	opacity =  DUSTTOGAS*100.*kappa(temper[l],dens[l]);  //ORIGINAL
	opar[l] = opacity; //Rosseland mean opacity (for radiative cooling)
	opap[l] = opacity; //Planck mean opacity (for radiative energy
			   //and thermal energy), it is chosen to be
			   //the same as Rosseland because they do not
			   //differ much in the temperature region we
			   //are interested in. (See papers of Bitsch
			   //et al 2013, 2014)
	opas[l] = DUSTTOGAS*100.*3.5*RHO0*R0;  //0.1*opacity; this is the stellar opacity	
      }
    }
  }
}

void RadEnergyToTemperature () 
{
  long gncell[3], stride[3];
  long i,j,h,l;
  real *temper, *energyrad,*eta1,*eta2;
  getgridsize (CurrentFluidPatch->desc, gncell, stride);
  temper    = CurrentFluidPatch->Temperature;
  energyrad = CurrentFluidPatch->EnergyRad;
  eta1      =  CurrentFluidPatch->Eta1;
  eta2      =  CurrentFluidPatch->Eta2;
  for (h = Nghost[2]; h < gncell[2]-Nghost[2]; h++) {
    for (i = Nghost[1]; i < gncell[1]-Nghost[1]; i++) {
      for (j = Nghost[0]; j < gncell[0]-Nghost[0]; j++) {
	l = j+i*stride[1]+h*stride[2];
	temper[l] = eta1[l] + energyrad[l]*eta2[l];
      }
    }
  }
  ExecCommSameOneField (CurrentFluidPatch->desc->level, temper); 
}

void TemperatureToEnergy () 
{
  long gncell[3], stride[3];
  long i,j,h,l;
  real *temper, *energy, *dens;
  getgridsize (CurrentFluidPatch->desc, gncell, stride);
  temper    = CurrentFluidPatch->Temperature;
  energy    = CurrentFluidPatch->Energy;
  dens      = CurrentFluidPatch->Density;
  for (h = 0; h < gncell[2]; h++) {
    for (i = 0; i < gncell[1]; i++) {
      for (j = 0; j < gncell[0]; j++) {
	l = j+i*stride[1]+h*stride[2];
#ifdef TEMPFLOOR
	if (temper[l] > BOUNDTEMP/TEMP0/TEMPRATIO)
#endif
	  energy[l] = CV*dens[l]*temper[l];
#ifdef TEMPFLOOR
	else
	  energy[l] = CV*dens[l]*BOUNDTEMP/TEMP0/TEMPRATIO;
#endif
      }
    }
  }
}

void ComputeRadiativeEnergy(real dt)
{
  FluidWork *fw;
  long size[3], gncell[3], stride[3];
  long i,j,h,l,nr,ns,ni;
  real *dens, *stellarrad, *opap, *colat, *rad;
  real *temper, *energyrad, *qplus,*eta1,*eta2;
  real denom, highcol=1e30;
  boolean TrueBoundary;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);

  temper    = fw->Temperature;
  energyrad = fw->EnergyRad;
  colat     = fw->desc->Center[_COLAT_];
  rad       = fw->desc->Center[_RAD_];

  if (!HalfDisk)
    highcol = RANGE3HIGH;
  
  for (i = 0; i < 3; i++)
    size[i] = fw->desc->gncell[i];
  
  nr = size[1]; // if coordpermut=213 (first coord is azimuth, 2nd rad, 3rd colat)
  ns = size[0]; // azimuth
  ni = size[2]; // co-lat
  
  dens = fw->Density;
  stellarrad =  fw->StellarRad;
  opap =  fw->OpaP;
  qplus=  fw->Qplus;
  eta1 =  fw->Eta1;
  eta2 =  fw->Eta2;
  

  // Eta is a quantity which connects the radiation energy to the
  // temperature. It depends on the viscous heating and the stellar
  // irradiation too. They are described in Bitsch et al. 2013 A&A 549,
  // 124 at equation B.7, B.8
  // (http://www.aanda.org/articles/aa/pdf/2013/01/aa20159-12.pdf)
  if (Restart == YES) FreshStart = NO;
  for (h = 0; h < ni; h++) {
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	TrueBoundary = NO;
	l = j+i*stride[1]+h*stride[2];
        denom = 1.0+16.0*dt*opap[l]*SIGMARAD*pow(temper[l],3.0)/CV;
	eta1[l] = temper[l]+12.0*dt*opap[l]/CV*SIGMARAD*pow(temper[l],4.0)+ \
	  dt*(stellarrad[l]+qplus[l])/(dens[l]*CV);
	eta1[l] = eta1[l]/denom;
	eta2[l] = dt*C0*opap[l]/(CV*denom);	
	
	if (rad[l] < RANGE2LOW) TrueBoundary = YES;
	if (rad[l] > RANGE2HIGH) TrueBoundary = YES;
	if (colat[l] < RANGE3LOW) TrueBoundary = YES;
	if (colat[l] > RANGE3HIGH) TrueBoundary = YES;
	if ((FreshStart == YES) || (TrueBoundary == YES))
	  energyrad[l] = (temper[l] - eta1[l])/eta2[l];
	if (colat[l] > RANGE3HIGH)
	  energyrad[l] = energyrad[l-stride[2]];
	if ((colat[l] < RANGE3LOW) || (colat[l] > highcol))
	  energyrad[l] = BoundEcode;
      }
    }
  }
  ExecCommSameOneField (fw->desc->level, energyrad); 
  ExecCommSameOneField (fw->desc->level, eta1); 
  ExecCommSameOneField (fw->desc->level, eta2); 
}
