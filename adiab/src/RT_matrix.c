#include "jupiter.h"

void GetRadEnergyDifference (real dt) {
  FluidWork *fw;
  long i,j,h,l,gncell[3],stride[3];
  real *diffenrad,*enrad,*enradold;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  enrad = fw->EnergyRad;
  enradold = fw->EnradOld;
  diffenrad = fw->DiffEnrad;
  for ( h = 0; h < gncell[2]; h++) { 
    for (i = 0; i < gncell[1]; i++) {
      for (j = 0; j < gncell[0]; j++) {
	l = j+i*stride[1]+h*stride[2];
	diffenrad[l] = (enrad[l]-enradold[l])/dt;
      }
    }
  }
  fw->Fluid->PreviousEradExists = YES;
}

void PredictNewRadEnergy (real dt) {
  FluidWork *fw;
  long i,j,h,l,gncell[3],stride[3];
  real *diffenrad,*enrad,*enradold;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  enrad = fw->EnergyRad;
  enradold = fw->EnradOld;
  diffenrad = fw->DiffEnrad;
  if (fw->Fluid->PreviousEradExists == YES) {
    for (h = 0; h< gncell[2]; h++) {
      for (i = 0; i < gncell[1]; i++) {
	for (j = 0; j < gncell[0]; j++){
	  l = j+i*stride[1]+h*stride[2];
	  enradold[l] = enrad[l];
	  enrad[l] += diffenrad[l]*dt;
	}
      }
    }
  } else {
    memcpy (enradold, enrad, gncell[0]*gncell[1]*gncell[2]*sizeof(real));
  }
}

void SolveMatrix(boolean everywhere)
{
  FluidWork *fw;
  long gncell[3], stride[3];
  long *pc;
  long i,j,h,l,lip,lim,ljp,ljm,lhp,lhm;
  int count,parity;
  real *rhs,*enrad;
  real *aar,*aap,*aat,*ccr,*cct,*ccp,*bb;
  char *hidden;
  real omega=1.0,anormrhs=0.0,anormrhstot=0.0,anorm=0.0,anorm0=1.0,anormtot=0.0,resid;
  real drad,dphi,dtheta,rhojacobi;
  real *center[3];
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  for (i = 0; i < 3; i++)	/* 3, not NDIM */
    center[i] = fw->desc->Center[i];
  pc = fw->desc->pcorner_min;


  enrad = fw->EnergyRad;
  aar = fw->Aarmatrix;
  aat = fw->Aatmatrix;
  aap = fw->Aapmatrix;
  ccr = fw->Ccrmatrix;
  cct = fw->Cctmatrix;
  ccp = fw->Ccpmatrix;
  bb  = fw->Bbmatrix;
  rhs = fw->Rhsmatrix;

  hidden = fw->desc->Hidden;

// Calculating the Norm of the Right-hand-side:
  anormrhs = 0.;
  for( h = Nghost[2]; h < gncell[2]-Nghost[2]; h++) {
    for (i = Nghost[1]; i < gncell[1]-Nghost[1]; i++) {
      for (j = Nghost[0]; j < gncell[0]-Nghost[0]; j++) {
    	l = j+i*stride[1]+h*stride[2];
	anormrhs += fabs(rhs[l]);
      }
    }
  }


// calculate spectral radius (rhojacobi, see Numerical recipies, page 867, Eq 19.5.24 spectral radius), THERE IS AN ERROR HERE IN PARALLEL, IT SHOULD BE NR, NS, NI FOR THE WHOLE SIMULATION NOT ONLY ON THE GIVEN CPU-PATCH, SO THESE NUMBERS SHOULD STAND FOR THE GLOBAL CELL NUMBERS NOT THE CPU PATCH CELL NUMBERS:

// Something better is needed

  	drad=fw->desc->Edges[_RAD_][1]-fw->desc->Edges[_RAD_][0];
	dphi= (fw->desc->Edges[_COLAT_][1]-fw->desc->Edges[_COLAT_][0])*RANGE2LOW;
	dtheta=(fw->desc->Edges[_AZIM_][1]-fw->desc->Edges[_AZIM_][0])*RANGE2LOW;
	rhojacobi=cos(M_PI/(real)SIZE2)+pow((drad/dphi),2.0)*cos(M_PI/(real)SIZE3)+pow((drad/dtheta),2.0)*cos(2*M_PI/(real)SIZE1);
	rhojacobi /= (1+pow((drad/dphi),2.0)+pow((drad/dtheta),2.0));


 // Compute the norm of rhs and distribute the result back to all processes: 
  MPI_Allreduce (&anormrhs, &anormrhstot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  resid = 0.;
  anormtot = anormrhstot;

  count = 0;
  while(anormtot > EPSILON*anormrhstot && count <20000){
    count ++;
    anorm = 0.0;
    
    // IMPORTANT: the limits of the following loops should go Ngh-1 till
    // gncell-1, otherwise the result is very different and wrong. If one
    // use only the active cells, the cooling of the disk stops quite
    // quickly, and do not evolve anymore... We should still understand
    // why...
    for( parity = 0; parity < 2; parity++) {
      for( h = Nghost[2]; h < gncell[2]-Nghost[2]; h++) {    
	for( i = Nghost[1]; i < gncell[1]-Nghost[1]; i++) {
	  for (j = Nghost[0]; j < gncell[0]-Nghost[0]; j++) {
	    if ((i+j+h+pc[0]+pc[1]+pc[2]) % 2 == parity) {
	      l = j+i*stride[1]+h*stride[2];
	      if ((hidden[l] != YES) || (everywhere == YES)) {
		lim = l-stride[1];
		lip = l+stride[1];
		ljp = l+1;
		ljm = l-1;
		lhp = l+stride[2];
		lhm = l-stride[2];
		resid = aar[l]*enrad[lim]+ccr[l]*enrad[lip]+	\
		  aat[l]*enrad[ljm]+cct[l]*enrad[ljp]+		\
		  bb[l]*enrad[l]+aap[l]*enrad[lhm]+		\
		  ccp[l]*enrad[lhp]-rhs[l];
		anorm += fabs(resid);
		enrad[l] -= omega*resid/bb[l];
	      }
	    }
	  }
	}
      }
      ExecCommSameOneField (fw->desc->level, enrad); 
    }
    if (count == 1) anorm0 = anorm;
    
    // Communicating the radiation energy Note: the comm below is needed
    // also for one CPU if the mesh is periodic in azimuth (as it
    // usually is)....
    
    if (count==1) 
      omega=1.0/(1.0-0.5*rhojacobi*rhojacobi);
    else 
      omega=1.0/(1.0-0.25*rhojacobi*rhojacobi*omega);

    MPI_Allreduce(&anorm, &anormtot, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    // End of iterative loop:
  }
  pInfo ("RT SOR counts: %d\n", count);

// Just a check printf:
/*  if (count < 0) {
    printf("Inside %s:\n", __FILE__);
    INSPECT_REAL (omega);
    INSPECT_REAL (anorm);
    INSPECT_REAL (anormtot);
    INSPECT_REAL (anormrhs);
    INSPECT_REAL (anormrhstot);
    INSPECT_INT  (count);
  }*/
}



