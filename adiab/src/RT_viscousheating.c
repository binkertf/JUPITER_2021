#include "jupiter.h"

void ComputeQplus ()
{
  FluidWork *fw;
  long i,j,h,l,lim,lip,ljp,ljm,lhp,lhm,gncell[3],stride[3];
  real heat=0.0,divergence;
  real *rho,*qplus,*rad,*theta,*phi,*vrad,*vtheta,*vphi;
  real Viscosity=0.0;
  fw = CurrentFluidPatch;
  if (fw->desc->level < VISCUTOFFLEVEL)
    Viscosity = VISCOSITY;
  getgridsize (fw->desc, gncell, stride);
  
  rho    = fw->Density;
  qplus  = fw->Qplus;
  rad    = fw->desc->Center[_RAD_];
  theta  = fw->desc->Center[_COLAT_];
  phi    = fw->desc->Center[_AZIM_];
  vrad   = fw->Velocity[_RAD_];
  vtheta = fw->Velocity[_COLAT_];
  vphi   = fw->Velocity[_AZIM_];

  for (l = 0; l < gncell[0]*gncell[1]*gncell[2]; l++) qplus[l] = 0.0;

  // ExecCommSameOneField (fw->desc->level, vrad);
  //  ExecCommSameOneField (fw->desc->level, vtheta);
  //  ExecCommSameOneField (fw->desc->level, vphi);
  //  ExecCommSameOneField (fw->desc->level, rho);

  //Please remember that : here theta is colatitude !!
  // and : vtheta and vphi are *angular*, not linear
  for (h = 1; h < gncell[2]-1; h++) {
    for (i = 1; i < gncell[1]-1; i++) {
      for (j = 1; j < gncell[0]-1; j++) {
	l = j+i*stride[1]+h*stride[2];
        lim = l-stride[1];
	lip = l+stride[1];
	ljp = l+1;
	ljm = l-1;
	lhp = l+stride[2];
	lhm = l-stride[2];
	divergence = (rad[lip]*rad[lip]*vrad[lip]-rad[lim]*rad[lim]*vrad[lim])/(rad[l]*rad[l]*	\
		     (rad[lip]-rad[lim]));
	divergence += (vtheta[lhp]*sin(theta[lhp])-vtheta[lhm]*sin(theta[lhm]))/      \
		      (sin(theta[l])*(theta[lhp]-theta[lhm]));
	divergence += (vphi[ljp]-vphi[ljm])/(phi[ljp]-phi[ljm]);
	heat  = 2.0*pow((vrad[lip]-vrad[lim])/(rad[lip]-rad[lim]),2.0);
	heat += 2.0*pow((vtheta[lhp]-vtheta[lhm])/(theta[lhp]-theta[lhm])+vrad[l]/rad[l],2.0);
	heat += 2.0*pow((vphi[ljp]-vphi[ljm])/(phi[ljp]-phi[ljm])+	\
	  		vrad[l]/rad[l]+vtheta[l]/tan(theta[l]),2.0);
	heat -= 2.0/3.0*pow(divergence,2.0);
      	heat += pow((rad[l]*(vtheta[lip]-vtheta[lim])/(rad[lip]-rad[lim]))+ \
      		    (vrad[lhp]-vrad[lhm])/(rad[l]*(theta[lhp]-theta[lhm])),2.0);
    	heat += pow(sin(theta[l])*(vphi[lhp]-vphi[lhm])/(theta[lhp]-theta[lhm])+(vtheta[ljp]-vtheta[ljm])/(sin(theta[l])*(phi[ljp]-phi[ljm])),2.0);
	heat +=  pow((vrad[ljp]-vrad[ljm])/(rad[l]*sin(theta[l])*(phi[ljp]-phi[ljm]))+rad[l]*sin(theta[l])* \
		     (vphi[lip]-vphi[lim])/(rad[lip]-rad[lim]),2.0);
	qplus[l] = heat*rho[l]*Viscosity;
        qplus[l] = 0.0; //THere is certainly a more clever way to get rid of this manual heating routine...
      }
    }
  }
  ExecCommSameOneField (fw->desc->level, qplus);

//WriteWorkArrayAlt(qplus, "qplus", 122);

}
