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
  int iwr;
  long gncell[3], i, j, k, m, stride[3];
  real *center[3]; 
  real SMOOTHORIG;
  real cellradius, limit, jup_2rad, const_dec;
  real xp,yp,zp,work1,work2,work3;
  real x,y,z,radius,azimuth,colatitude; /* have to  be used coherently
					   together   with  the  input
					   file  (e.g. do not  use 'x'
					   in  spherical  coordinates,
					   except   as   a  temporary,
					   unitialized variable */
  real pot;
  real mass = 0.0, tmass=0.0, cut_dens, cellvolume;
  FILE * sFile;
  getgridsize (fp->desc, gncell, stride);
  for (i = 0; i < 3; i++)	/* 3, not NDIM */
    center[i] = fp->desc->Center[i];

// ========================== modifications made from here: ============================

  SMOOTHING = (center[_RAD_][stride[1]]-center[_RAD_][0])*(center[_AZIM_][1]-center[_AZIM_][0])*(center[_COLAT_][stride[2]]-center[_COLAT_][0]);
  cellvolume= (center[_RAD_][stride[1]]-center[_RAD_][0])*(center[_AZIM_][1]-center[_AZIM_][0])*(center[_COLAT_][stride[2]]-center[_COLAT_][0]);
  SMOOTHING = 2.*pow(SMOOTHING,1./3.);
  SMOOTHORIG=SMOOTHING;

  jup_2rad=1.7959e-04;
 
  if(SMOOTHING < jup_2rad){ 
    limit=jup_2rad;
  } else {
    limit=SMOOTHING;
  }
  
  if(GlobalDate<=GlobalDateInit){
//    printf("at Globaldate %.20g level %ld smoothing %lg\n",GlobalDate,fp->desc->level,limit);
    iwr=0;
  }

  if (fp->desc->level == FinestLevel){
    if(GlobalDate<2*3.14159+GlobalDateInit){
      SMOOTHING = (2.*SMOOTHORIG-limit) * pow(cos((GlobalDate-GlobalDateInit)/4.),2)+limit;
//SMOOTHING = (2.*SMOOTHORIG-limit) * pow(cos((GlobalDate-GlobalDateInit)/24.),2)+limit;
//    if(sin(GlobalDate-GlobalDateInit)>=0.){
//      if(iwr==0)printf("in level 6 at date %.20g smothing %lg\n",GlobalDate,SMOOTHING);
//      iwr=1;}else{iwr=0;}
    }
  }

  if(SMOOTHING < limit) SMOOTHING=limit;
 /* if (fp->desc->level == 3) SMOOTHING*=2.0;
  if (fp->desc->level == 4) SMOOTHING*=2.0;
  if (fp->desc->level == 5) SMOOTHING*=2.0;
  if (fp->desc->level == 6) SMOOTHING*=2.0;*/

    pInfo ("Smoothing at %lg on lev %d: %lg\n", GlobalDate, fp->desc->level, SMOOTHING);


// ============================ till here ============================

  for (k = 0; k < gncell[2]; k++) { /* ghosts are included */
    for (j = 0; j < gncell[1]; j++) {
      for (i = 0; i < gncell[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
// modification made from here
/*        if  (CurrentFluidPatch->Density[m]> 1.422693105581794E+05) {
	  if (fp->desc->level == FinestLevel && k > 1 && k < gncell[2]-2 && j >1 && j <  gncell[1]-2 && i >1 && i <  gncell[0]-2){
	    cut_dens=CurrentFluidPatch->Density[m]-1.422693105581794E+05;
	    mass=mass+cut_dens*cellvolume;
	  }
	  CurrentFluidPatch->Density[m]=1.422693105581794E+05 ; 
        }*/
// till here

//       if  (CurrentFluidPatch->Density[m]>5.592728535953456E+05) {
//CurrentFluidPatch->Density[m]=5.592728535953456E+05 ;	
//}
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
  // modification for the sink-mass computation:
/*  MPI_Allreduce (&mass, &tmass, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if  (CPU_Rank == 0 && fp->desc->level == FinestLevel && tmass > 0.0 ) {
    sFile=fopen("/scratch/jszulagyi/1J_nodamp/accr_mass.dat","a+");
    fprintf(sFile,"%.10g\t%.15g\n",tmass,GlobalDate);
    fclose(sFile);
  }*/
// till here
}
