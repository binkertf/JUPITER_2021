#include "jupiter.h"

void SetRTBoundaryConditions () {
  FluidWork *fw;
  real *aar,*ccr,*ccp,*bb;
  real *rad,*col;
  long i,j,h,l,gncell[3],stride[3];
  real steprad,stepcol;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  rad = fw->desc->Center[_RAD_];
  steprad = fw->desc->Edges[_RAD_][1]-fw->desc->Edges[_RAD_][0];
  col = fw->desc->Center[_COLAT_];
  stepcol = fw->desc->Edges[_COLAT_][1]-fw->desc->Edges[_COLAT_][0];
  aar = fw->Aarmatrix;
  ccr = fw->Ccrmatrix;
  ccp = fw->Ccpmatrix;
  bb  = fw->Bbmatrix;
  //  Boundary condition for matrix	
  /* In case of Reflecting Boudaries */
  //We put a mirror for photons at the inner radius (edge of active mesh)
  for( h = 1; h < gncell[2]-1; h++) {
    for (i = 0; i <= Nghost[1]; i++) {
      for (j = 1; j < gncell[0]-1; j++) {
	      l = j+i*stride[1]+h*stride[2];
	      if(rad[l] < RANGE2LOW+steprad){
	        bb[l] +=  aar[l];
	        aar[l] = 0.;
	      }
      }
    }
  }
  
  //We put a mirror for photons at the outer radius (edge of active mesh)
  for( h = 1; h < gncell[2]-1; h++) {
    for (i = gncell[1]-1; i >= gncell[1]-Nghost[1]-1; i--) {
      for (j = 1; j < gncell[0]-1; j++) {
	      l = j+i*stride[1]+h*stride[2];
	      if(rad[l] > RANGE2HIGH-steprad){
	        bb[l] +=  ccr[l];
	        ccr[l] = 0.;
	      }
      }
    }
  }
  //If we simulate only half the disk, we put a mirror at the equator.
  //  if(HalfDisk){
    for (h = gncell[2]-1; h >= gncell[2]-Nghost[2]-1; h--) {
      for (i = 1; i < gncell[1]-1; i++) {
	      for (j = 1; j < gncell[0]-1; j++) {
	        l = j+i*stride[1]+h*stride[2];
	        if(col[l] > RANGE3HIGH-stepcol) {
	          bb[l] += ccp[l];
	          ccp[l] = 0;
	        }
	      }
      }
    }
  //}
}
