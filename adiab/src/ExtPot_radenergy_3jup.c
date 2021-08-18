/* 

$Author: masset $

$Id: ExtPot.c,v 5.24 2008/06/03 16:05:53 masset Exp $

$Log: ExtPot.c,v $
Revision 5.24  2008/06/03 16:05:53  masset
Multifluide works except coupling and flux communications


This file corresponds to  the implementation of an external potential,
that is to say an  arbitrary potential with an analytic dependence (it
may  be  time dependent),  not  necessarily  solution  of the  Poisson
equation (the system is  not self-gravitating).  This potential can be
either  specified in  the potential  library,  which is  found in  the
'libpot.txt' file. In this case, the parameter file should include the
parameter 'PotentialCode', set to the string identifying the potential
set up  (e.g. 'kepler').  Note that  you can add any  potential law to
this  library  (it  may  contain   up  to  511  different  laws,  edit
includes/potcodes.h to change this maximum value). The syntax is quite
simple, one should simply be  careful not to forget the semicolumn (;)
at the  end of  each line, except  those containing the  potential law
string identifier,  which end with a  columnb (:). At  build time, the
libpot.txt file  is processed by the perl  script 'potparser.pl' which
creates  automatically 'libpot.cx', 'potcode.c',  and 'includes/pot.h'
(and therefore overrides any  manual change to these files). Basically
this  script  does  the  following  :  (1)  it  converts  the  library
libpot.txt  into a  C source,  in the  file libpot.cx;  note  the 'cx'
suffix  instead of  'c'.   This is  meant  to ensure  that  it is  not
included in the  CVS distribution through the use  of wildcards (*.c),
since  this  file can  be  built  using the  files  given  by the  CVS
distribution.   (2)  It  converts  the  law codes  into  integers  (in
includes/pot.h) and  (3) it associates  the integer to each  string in
potcode.c. Note  that the code will  go through this  procedure if and
only if the  boolean parameter 'ExternalPotential' is set  to YES (see
'actual.c'). It defaults  to NO (see 'var.c'). It  is nevertheless set
automatically to YES if  'PotentialCode' is defined (see switch.c), so
there is  no need to  specify explicitely 'ExternalPotential'  in this
case. Note  that you can also  specify manually the  expression of the
potential  in the present  source file,  between the  comments clearly
stating where  your definition  should be included.  In this  case you
want to remove any potential  code from your parameter file, otherwise
the code will detect that the potential law is multiply defined and it
will issue an error. */

/* Note for writing  a potential law in libpot.txt:  all the variables
x,y,z,radius,azimuth  and colatitude  are defined,  regardless  of the
mesh  geometry.  This  enables  one to  use  e.g. 'x'  as a  temporary
variable say  in cylindrical or spherical geometry,  without having to
declare another variable below */

#include "jupiter.h"
#include "pot.h"

static boolean NeverTested = TRUE;

void ComputeExternalPotential (GlobalDate, fp, t, phi, flag)
     FluidPatch *fp;
     real t;
     real *phi;
     char flag;
     real GlobalDate;
{
  long gncell[3], i, j, k, m, stride[3];
  real *center[3];
  real SMOOTHORIG;
  real limit, jup_2rad;
  real  cut_dens,mass=0.0,cellvolume, tmass=0.0,cut_dens2,mass2=0.0, tmass2=0.0;
  real temp1=0.0,ttemp1=0.0,temp2=0.0,ttemp2=0.0;
  real dr,daz,dco;
  FILE* sFile;
  real *dens, *erad,*temp, *InvVolume;
  real erad_hill=0.0,erad_lev6=0.0, erad_lev3=0.0,terad_hill=0.0,terad_lev6=0.0, terad_lev3=0.0;
  real xp=1.0,yp=0.0,zp=0.0,work1,work2,xr,yr,zr,dist,distx,disty,distz; 
  real x,y,z,radius,azimuth,colatitude; /* have to  be used coherently
					   together   with  the  input
					   file  (e.g. do not  use 'x'
					   in  spherical  coordinates,
					   except   as   a  temporary,
					   unitialized variable */
  real pot;
  dens = fp->Density->Field;  // desnity is called for the densityfloor
  erad = fp->EnergyRad->Field;  // radiative energy
  temp = fp->Temperature->Field; 
  InvVolume = fp->desc->InvVolume;
  getgridsize (fp->desc, gncell, stride);
  for (i = 0; i < 3; i++)	/* 3, not NDIM */
    center[i] = fp->desc->Center[i];

// ========================== modifications made from here: ============================
// Smoothing-tapering:



  //  SMOOTHING =  sqrt((center[_RAD_][stride[1]]-center[_RAD_][0])*(center[_RAD_][stride[1]]-center[_RAD_][0])+\
		(center[_AZIM_][1]-center[_AZIM_][0])*(center[_AZIM_][1]-center[_AZIM_][0])+ \
      		(center[_COLAT_][stride[2]]-center[_COLAT_][0])*(center[_COLAT_][stride[2]]-center[_COLAT_][0]));
    cellvolume= (center[_RAD_][stride[1]]-center[_RAD_][0])*(center[_AZIM_][1]-center[_AZIM_][0])* \
      (center[_COLAT_][stride[2]]-center[_COLAT_][0]);
 //   SMOOTHING = SMRATIO*SMOOTHING; // diagonal of a cell!
 /*   SMOOTHING = (center[_RAD_][stride[1]]-center[_RAD_][0])*(center[_AZIM_][1]-center[_AZIM_][0])* \
      (center[_COLAT_][stride[2]]-center[_COLAT_][0]); // old fashion smoothing
    SMOOTHING = SMRATIO*pow(SMOOTHING,1./3.); // diameter of TWO cell!*/

// Elena's smoothing lengths:
 if (fp->desc->level == 0)  SMOOTHING=8*0.003379068;
 if (fp->desc->level == 1)  SMOOTHING=4*0.003379068;
 if (fp->desc->level == 2)  SMOOTHING=2*0.003379068;
 if (fp->desc->level == 3)  SMOOTHING=0.003379068;
 if (fp->desc->level == 4)  SMOOTHING=0.003379068;
 if (fp->desc->level == 5)  SMOOTHING=0.5*0.003379068;
// if (fp->desc->level == 4)  SMOOTHING=0.5*0.003379068;
// if (fp->desc->level == 5)  SMOOTHING=0.25*0.003379068;
// if (fp->desc->level == 6)  SMOOTHING=0.125*0.003379068;
 if (fp->desc->level == 6)  SMOOTHING=0.5*0.003379068;

    SMOOTHORIG=SMOOTHING;
    
    
    jup_2rad=1.7959e-04;  /* this is the diameter of Jupiter (the
			  planet's...), if the smoothing-length would
			  be smaller than this, we use the planet's
			  diameter as the smoothing length */
    if(SMOOTHING < jup_2rad){ 
      limit=jup_2rad;
    } else {
      limit=SMOOTHING;
    }
    
    
    if (SmoothTaper == YES && fp->desc->level == LevMax){
      if(GlobalDate<2*3.14159+GlobalDateInit){
	// smooth-taper for half orbit of the planet
	SMOOTHING = (2.*SMOOTHORIG-limit) * pow(cos((GlobalDate-GlobalDateInit)/4.),2)+limit; 
      }
    }
    if(SMOOTHING < limit) SMOOTHING=limit;

 //   if (fp->desc->level == 1)  SMOOTHING=SMOOTHING*2.;
   // if (fp->desc->level == 2)  SMOOTHING=SMOOTHING*2.;

//    if (fp->desc->level ==5)  SMOOTHING=SMOOTHING*2.;
  //  if (fp->desc->level ==6)  SMOOTHING=SMOOTHING*4.;
    pInfo ("Smoothing at %g on lev %d: %g\n", GlobalDate, fp->desc->level, SMOOTHING);


// ============================ till here ============================
  for (k = 0; k < gncell[2]; k++) { /* ghosts are included */
    for (j = 0; j < gncell[1]; j++) {
      for (i = 0; i < gncell[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];

	if (fp->desc->CommSource[m] & flag) {
	  radius     = x  = center[_RAD_][m];
	  azimuth    = y  = center[_AZIM_][m];
	  colatitude = z  = center[_COLAT_][m];
	  pot = -1e29;

	  /* User needs to edit the following line(s) for his own use */
	  
	  /****************************************************************/
	  /*                                                              */
	  /*       Please Edit BELOW and only BELOW this comment          */
	  /*                                                              */
	  /****************************************************************/
	  /* If you want to paste a patch please do it below this comment */
	  /* and remove the present potential set up                      */
	  /****************************************************************/


	  /****************************************************************/
	  /* Use any of the following: x,y,z,radius,azimuth,colatitude,t  */
	  /****************************************************************/
	  /*                                                              */
	  /*       Please Edit ABOVE and only ABOVE this comment          */
	  /*                                                              */
	  /****************************************************************/
	  
	  /* After this line the source file should NOT be edited */
#include "libpot.cx"
	  phi[m] = pot;
	}
      }
    }
  }
// modification made from here
// This part is for monitoring the mass in the innermost cells:
   for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];

        if  (fp->desc->level == LevMax){
	    dr=center[_RAD_][stride[1]]-center[_RAD_][0];
	    daz=center[_AZIM_][1]-center[_AZIM_][0];
	    dco=center[_COLAT_][stride[2]]-center[_COLAT_][0];
	    if(center[_AZIM_][m] < 0.0+daz*1.2 && center[_AZIM_][m] > 0.0-(daz*1.2) && center[_RAD_][m] < 1.0+dr*1.2 &&  center[_RAD_][m] > 1.0-(dr*1.2) && center[_COLAT_][m] > (M_PI/2.)-(dco*1.2) &&  center[_COLAT_][m] < (M_PI/2.)) {
//cut down density to enforce accretion:
/*     		  if  (CurrentFluidPatch->Density[m]> 4.e+01) {
	   	  cut_dens=CurrentFluidPatch->Density[m]-4.e+01;
	    	  mass=mass+cut_dens/InvVolume[m];
	  	  CurrentFluidPatch->Density[m]=4.e+01 ; 
       	          }*/
	    cut_dens=CurrentFluidPatch->Density[m];
	    mass=mass+cut_dens/InvVolume[m];
	    temp1=temp1+temp[m];
	    }
    	    if(center[_AZIM_][m] < 0.0+daz*2.2 && center[_AZIM_][m] > 0.0-(daz*2.2) && center[_RAD_][m] < 1.0+dr*2.2 &&  center[_RAD_][m] > 1.0-(dr*2.2) && center[_COLAT_][m] > (M_PI/2.)-(dco*2.2) &&  center[_COLAT_][m] < (M_PI/2.)) {
	    cut_dens2=CurrentFluidPatch->Density[m];
	    mass2=mass2+cut_dens2/InvVolume[m];
	    temp2=temp2+temp[m];
	    }
        }
// Monitor the radiative energy fluctuations during the kappa mechanism:
			xr = center[_RAD_][m]*cos(center[_AZIM_][m])*sin(center[_COLAT_][m]);
			yr = center[_RAD_][m]*sin(center[_AZIM_][m])*sin(center[_COLAT_][m]);
			zr = center[_RAD_][m]*cos(center[_COLAT_][m]);
			distx=(xr-xp);
			disty=(yr-yp);
			distz=(zr-zp);
			dist=sqrt(distx*distx+disty*disty+distz*distz);	
	if(dist < 0.1) {
	    erad_hill=erad_hill+erad[m];
	}
	if(fp->desc->level == 3) {
	    erad_lev3=erad_lev3+erad[m];
	}
	if(fp->desc->level == 6) {
	    erad_lev6=erad_lev6+erad[m];
	}
      }
    }
  }
  MPI_Allreduce (&mass, &tmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&mass2, &tmass2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&temp1, &ttemp1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&temp2, &ttemp2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);  
  MPI_Allreduce (&erad_hill, &terad_hill, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&erad_lev3, &terad_lev3, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce (&erad_lev6, &terad_lev6, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if  (CPU_Rank == 0 && fp->desc->level == LevMax ) {
  pInfo ("Planet mass %g, mass with 2 cell radius: %g at date %.15g\n", 2.0*tmass,2.0*tmass2, GlobalDate);
  pInfo ("Temperature average 4 cells %g, 32 cells: %g at date %.15g\n", ttemp1/4.0,ttemp2/32., GlobalDate);
  pInfo ("Erad integral Hill-sphere %g, Level 3: %g, Level 6: %g at date %.15g\n",terad_hill*2.0, erad_lev3*2.0,erad_lev6*2.0, GlobalDate);
  }
// till here
}
