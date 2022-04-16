#include "jupiter.h"
#include "calc_gamma.h"

#define MAXLEVELIRRAD 100
static MPI_Comm RadialCPUBeam[MAXLEVELIRRAD];
static MPI_Comm HorizontalCPUSlice[MAXLEVELIRRAD];
static boolean  CommCreated[MAXLEVELIRRAD];

static real InnerRings_local[2048];
static real InnerRings_reduced[2048];

void  ComputeStellarHeating()
{
  FluidWork *fw;
  long gncell[3], stride[3];
  long i, j, h, l;
  real *tau,tauprev,*dens,*stellarrad,*opas;
  real *radint, *col, *rad;
  real fstarvol, deltar;

  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);

  dens       = fw->Density;
  tau        = fw->TauOptical;
  stellarrad = fw->StellarRad;
  opas       = fw->OpaS;
  radint     = fw->desc->Edges[_RAD_];
  rad        = fw->desc->Center[_RAD_];
  col        = fw->desc->Center[_COLAT_];
  deltar     = radint[1]-radint[0];

  ComputeOpticalDepth(); // here we call the computation of optical depth (tau)

  ExecCommSameOneField (fw->desc->level, tau);

  for (l = 0; l < gncell[0]*gncell[1]*gncell[2]; l++) stellarrad[l] = 0.0;

  //Calculating the stellar volume and surface and then calculating the
  for( h = 0; h < gncell[2]-Nghost[2]; h++) {
    if (col[h*stride[2]] > RANGE3LOW) {
      for( i = 0; i < gncell[1]; i++) {
	fstarvol = FSTAR/(radint[i]*radint[i]+radint[i]*deltar+deltar*deltar/3.0);
	for (j = 0; j < gncell[0]; j++) {
	  l = j+i*stride[1]+h*stride[2];
	  if (rad[l] > RANGE2LOW) { //We do not heat in the innermost ghosts
	    tauprev = tau[l]-dens[l]*opas[l]*(radint[i+1]-radint[i]); //We cannot use lim as it would not work for leftmost ghost.
	    if((col[l]/M_PI*180.0) <= (90.-STELLDEG) || (col[l]/M_PI*180.0) >= (90.+STELLDEG)) {
	      //only stellar irradiation for disc above the shadow
	      stellarrad[l] = fstarvol*exp(-tauprev)*(1.0-exp(-dens[l]*opas[l]*deltar))/(deltar);
	    }
	    else {
	      stellarrad[l] = 0.0;
	    }
	  }
	}
      }
    }
  }
 ExecCommSameOneField (fw->desc->level, stellarrad);
}



void ComputeGammaField() {
  
  long i,j,h,l,nr,ns,ni;
  real *dens, *temp,*gamma_out;
  long stride[3],gncell[3];
  FluidWork *fw;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);

  nr = gncell[1];
  ns = gncell[0];
  ni = gncell[2];

  dens = fw->Density;
  temp = fw->Temperature;
  gamma_out = fw->Gamma;

  for( h = 0; h < ni; h++) {
    for ( i = 0; i < nr; i++ ) {
      for ( j = 0; j < ns; j++ ) {
	      l = j+i*stride[1]+h*stride[2];
        gamma_out[l] = Gamma1(temp[l], dens[l]);
  
      }
    }
  }
  
}


void CreateRadialComm () {
  int level;
  tGrid_CPU *desc;
  desc = CurrentFluidPatch->desc;
  level = desc->level;
  if (!CommCreated[level]) {
    MPI_Comm_split (MPI_COMM_WORLD, desc->color, desc->key, &RadialCPUBeam[level]);
    MPI_Comm_split (MPI_COMM_WORLD, desc->colorz, 0, &HorizontalCPUSlice[level]);
    CommCreated[level] = YES;
  }
}

void ComputeOpticalDepth()
{
  FluidWork *fw;
  long size[3], gncell[3], stride[3];
  real *tau,*dens, *opas, *taumax, *depth, *radint, *rad, *tau_out;
  long i, j, h, l, nr, ns, ni, lim, l2D, lev, key, nrlast;
  int CountUs_0, CountUs;
  // CountUs: A counter to count the CPUs which touch
  // the inner boundary in a slice of given altitude
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  radint = fw->desc->Edges[_RAD_];
  rad    = fw->desc->Center[_RAD_];
  depth  = fw->desc->optical_depth;
  lev    = fw->desc->level;
  key    = fw->desc->key;

  for (i = 0; i < 3; i++)
    size[i] = fw->desc->gncell[i];

  nr = size[1]; // if coordpermut=213 (first coord is azimuth, 2nd rad, 3rd colat)
  ns = size[0]; // azimuth
  ni = size[2]; // co-lat

  nrlast = nr-Nghost[1]-1;

  tau = fw->TauOptical;
  taumax = fw->TauMax;
  dens = fw->Density;
  opas = fw->OpaS;
  tau_out = fw->TauCell;


  // Average the densities azimuthally to avoid perturbations:
 i=Nghost[1]-1;
  for(h =0; h< ni; h++) {
    InnerRings_local[h]=0.0;
    for (j = 0; j < ns; j++) {
      l = j+i*stride[1]+h*stride[2];
      if (rad[l] < RANGE2LOW)
	InnerRings_local[h] += dens[l];
    }
    InnerRings_local[h] /= (real)ns;
    if (rad[l] < RANGE2LOW) {
        CountUs=1;
    } else {
        CountUs=0;
    }
  }

  CreateRadialComm ();
  /* Watch out: the sums below must be done at a given altitude: if we
     have two rows of CPUs vertically, we cannot sum up the density
     zones of the lower row with those of the upper row: this would be
     meaningless (it will cut the radiation at high altitude since we
     add the high densities of the equator in the average...). This is
     why the communicator must not be MPI_COMM_WORLD but
     HorizontalCPUSlice[lev], which is a subset (level dependent) of
     the CPUs at same altitude (color: CPU altitude i[2] in
     split.c)  */
  MPI_Allreduce(InnerRings_local,InnerRings_reduced,ni,MPI_DOUBLE,MPI_SUM,HorizontalCPUSlice[lev]);
  MPI_Allreduce(&CountUs,&CountUs_0,1,MPI_INT,MPI_SUM,HorizontalCPUSlice[lev]);


  /* Note that the sentence below yields Nans for lev>0, since
   CountUs_0 = 0 for refined patches (they do not intersect the inner
   boundary */
  for(h =0; h < ni; h++) InnerRings_local[h]=InnerRings_reduced[h]/(real)CountUs_0;

  //////////////////////////////////////////////////////////////////////////
  // We define here the optical depth for the first ring (of the whole disk)
  // The radial inner ghost cells. This part is user dependent.
  i=Nghost[1]-1;  //It can be specified in the last line of ghosts
		  //only. Anything else is ignored.
  for(h =0; h< ni; h++) {
    for (j = 0; j < ns; j++) {
      l = j+i*stride[1]+h*stride[2];
      /* Note that refined patches never enter the following test, so
	 the Nans they have in InnerRings_local are harmless */
      if (rad[l] < RANGE2LOW)
	tau[l] = 2.0*InnerRings_local[h]*opas[l]*(radint[i+1]-radint[i]);
      // Using here the azimuthally averaged densities
    }
  }
  //End of definition of optical depth of ghost ring(s).
  //////////////////////////////////////////////////////////////////////////

  for (h = 0; h < ni; h++) {
    for (i = Nghost[1]; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l   = j+i*stride[1]+h*stride[2];
	l2D = j+h*gncell[0];
	lim = l-stride[1];
	tau[l] = dens[l]*opas[l]*(radint[i+1]-radint[i]); //Optical depth of the cell
	tau_out[l] = tau[l];
	if ((i > Nghost[1]) || (key == 0)) //We add the optical depth
					   //of the previous active
					   //cell (in r) if it exists,
					   //or that of the ghost if
					   //we are on the innermost
					   //processor
	  tau[l] += tau[lim];
	if (i == nrlast)
	  depth[l2D] = tau[l];
      }
    }
  }

  MPI_Scan (MPI_IN_PLACE, depth, gncell[0]*gncell[2], MPI_DOUBLE, MPI_SUM, RadialCPUBeam[lev]);

  //The sum above is done in place: it contains its own value on the
  //processor of (local) rank 0, the sum of ranks 0 and 1 on rank 1,
  //etc.

  for (h = 0; h < ni; h++) {
    for (i = 0; i < nr; i++) {
      for (j = 0; j < ns; j++) {
	l   = j+i*stride[1]+h*stride[2];
	l2D = j+h*gncell[0];
	tau[l] += depth[l2D]-tau[j+h*stride[2]+nrlast*stride[1]];
	//  We want to add the optical depth of the *previous* CPUs (in
	//  r). We therefore have to remove the sum 'depth' the optical
	//  depth of the current CPU (which has been taken into account
	//  in the MPI_Scan above). What we remove in the last term is
	//  exactly what depth[l2D] was set to prior to MPI_Scan().
      }
    }
  }
}
