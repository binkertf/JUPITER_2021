/*

$Author: masset $

$Id: Init.c,v 5.20 2008/06/10 15:51:18 masset Exp $

$Log: Init.c,v $
Revision 5.20  2008/06/10 15:51:18  masset
*** empty log message ***

Revision 5.19  2008/06/06 17:57:06  masset
Multifluid approximately OK. Scalability test (on rho) failed

Revision 5.18  2008/06/03 16:05:53  masset
Multifluide works except coupling and flux communications

Revision 5.17  2008/02/05 05:29:40  masset
High lev proj works but in seq built only

*/

#include "jupiter.h"
#include "init.h"

void InitialCondition (fp)
     FluidPatch *fp;
{
  real *energy_field, *density_field, *velocity[3], *center[3];
  long stride[3], gncell[3], dim, i, j, k, m;
  real x,y,z,radius,azimuth,colatitude; /* Will have to be used coherently together with the input file */
  real vx, vy, vz, vrad, v_azimuth, v_colatitude;
  real Bx, By, Bz, Brad, B_azimuth, B_colatitude;
  real density, energy, axial_radius;
  real interm1, interm2, interm3; /* Work variables for user convenience */
  /* They are unused in most setups */
  energy_field = fp->Energy->Field;
  density_field = fp->Density->Field;
  getgridsize (fp->desc, gncell, stride);
  for (dim = 0; dim < 3; dim++) { /* Do NOT change this 3 to NDIM (see below) */
    velocity[dim] = fp->Velocity->Field[dim];
    center[dim] = fp->desc->Center[dim];
  }
  for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
    for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
      for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
	m = i*stride[0]+j*stride[1]+k*stride[2];
	radius     = x  = center[_RAD_][m];
	azimuth    = y  = center[_AZIM_][m];
	colatitude = z  = center[_COLAT_][m];
	axial_radius = radius * sin(colatitude);

	density = energy = -1e20;
	axial_radius=0.0;
	vx = vy = vz = vrad = v_azimuth = v_colatitude = 0.0;
	Bx = By = Bz = Brad = B_azimuth = B_colatitude = 0.0;

	/* User needs only to edit the following lines for his own use */
	/* but note that it is a better practice to create an entry */
	/* in the file 'libinit.txt' */
	/* If you want to paste a patch please do it below this comment */
	/* and remove any prior initial condition set up that remains */
	/* between this comment and the next one */
	/***************************************************************/


	/***************************************************************/
	/* If you want to paste a patch please do it above this comment */
	/* and remove the lines corresponding to the present setup */
	/* Initialize ... (input: x,y,z,radius,azimuth,colatitude) */
	/* output: density, energy, vx, vy, vz, Bx, By, Bz etc.    */

	/* After this line the source file should NOT be edited */

#include "libinit.cx"

	switch (CoordType) {
	case CARTESIAN:
	  if (CoordNb[0] < NDIM) velocity[CoordNb[0]][m] += vx;
	  if (CoordNb[1] < NDIM) velocity[CoordNb[1]][m] += vy;
	  if (CoordNb[2] < NDIM) velocity[CoordNb[2]][m] += vz;
	  break;
	case CYLINDRICAL:
	  if (CoordNb[0] < NDIM) velocity[CoordNb[0]][m] += vrad;
	  if (CoordNb[1] < NDIM) velocity[CoordNb[1]][m] += v_azimuth;
	  if (CoordNb[2] < NDIM) velocity[CoordNb[2]][m] += vz;
	  break;
	case SPHERICAL:
	  if (CoordNb[0] < NDIM) velocity[CoordNb[0]][m] += vrad;
	  if (CoordNb[1] < NDIM) velocity[CoordNb[1]][m] += v_azimuth;
	  if (CoordNb[2] < NDIM) velocity[CoordNb[2]][m] += v_colatitude;
	  break;
	}
	density_field[m] += density;
	if (density_field[m] < _Small_Rho)
	  _Small_Rho = density_field[m];
	energy_field[m] += energy;
      }
    }
  }
}

void InitWholeHierarchy (NbRestart)
     long *NbRestart;
{
  tGrid_CPU *item;
  real small;
  long lev;
  static boolean DivideSmallRho = NO;
  FluidPatch *Fluid;
  real *density_field;
  long i,j,k,m,gncell[3],stride[3];
  item = Grid_CPU_list;
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      Fluid = item->Fluid;
      while (Fluid != NULL) {
	SetFluidProperties (Fluid);
//	if (IsothermalToStellar) {
//	   InitialCondition (Fluid);
//	   Stellar = NO;
//	   ReadField (Fluid, *NbRestart);
//	   Stellar = YES;
//	}
	if (!Restart) {
	  if (!Stretch)
	    InitialCondition (Fluid);
	  else {
	    if (Fluid->desc->level > 0) {
	      prs_error ("Array stretching allowed only on non-refined setups");
	    }
	    ReadField (Fluid, *NbRestart);
	    density_field = Fluid->Density->Field;
	    getgridsize (Fluid->desc, gncell, stride);
	    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
	      for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
		for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
		  m = i*stride[0]+j*stride[1]+k*stride[2];
		  if (density_field[m] < _Small_Rho)
		    _Small_Rho = density_field[m];
		}
	      }
	    }
	  }
	}
	else
	  	 ReadField(Fluid, *NbRestart);
	if (SuperImpose)
	  InitialCondition (Fluid);
	TrueBC_fp (Fluid);
	Fluid = Fluid->next;
      }
    }
    item = item->next;
  }
  MPI_Allreduce (&_Small_Rho, &small,  1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  if (!DivideSmallRho) {
    _Small_Rho = small*SMALLESTDENSRATIO;
    DivideSmallRho = YES;
  }
  _Small_Energy = _Small_Rho * ASPECTRATIO*ASPECTRATIO; // Not dimensionally homogeneous
  CommAll ();





  // The RT needs to be initialized from top to bottom.
  if (Stellar) { //We synchronize boundaries of Energy Rad fields prior to next loop. Needed for restart
    for (lev = 0; lev <= LevMax; lev++) {
      item = Grid_CPU_list;
      while (item != NULL) {
	if ((item->cpu == CPU_Rank) && (item->level == lev)) {
		Fluid = item->Fluid;
      	while (Fluid != NULL) {

	  SendToCurrent (Fluid);
		//printf("RT initalization: %s\n",Fluid->Name); 
	  ExecCommSameOneField (CurrentFluidPatch->desc->level, CurrentFluidPatch->EnergyRad);
	  CurrentToPatch (Fluid);
	Fluid = Fluid->next;

		}
	}
	item = item->next;
      }
    }
  }




  if (Stellar) {
    for (lev = 0; lev <= LevMax; lev++) {
      GetOpticalDepthFromLevel (lev-1);
      GetEnergyRadFromLevel (lev-1);
      item = Grid_CPU_list;
      while (item != NULL) {
	if ((item->cpu == CPU_Rank) && (item->level == lev)) {


		Fluid = item->Fluid;
      	while (Fluid != NULL) {


	  SendToCurrent (Fluid);
	  ComputeTemperatureField ();
	  ComputeOpacityInit ();
	  ComputeQplus ();
	  ComputeStellarHeating ();
	  if (Restart != YES)
	    ComputeRadiativeEnergy (DT/1.0e6);
	  ComputeDiffusionCoefficients ();
	  if (Restart != YES)
	    ComputeMatrixElements (DT/1.0e6);
	  SetRTBoundaryConditions ();
	  if (Restart != YES)
	    SolveMatrix (YES); // Solve elliptic eq everywhere (including above refined patches)
	  CurrentToPatch (Fluid);
		Fluid = Fluid->next;
		}
	}
	FluidCoupling (item, 0.0);
	item = item->next;
      }
    }
  }
  FreshStart = NO;
  if (Stretch == YES)
    *NbRestart = 0;
}
