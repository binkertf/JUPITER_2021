#include "jupiter.h"

void Predictor_adiab (dt)
     real dt;
{
  long i, j, k, l, m, n, p, Size;
  FluidWork *fw;
  real *rho, *rhop, *e, *ep, var, temperature;
  real *v[3], *vp[3], *source_v[3], *source_div;
  real *slope_rho[3], *slope_v[3][3], *slope_e[3];
  long gncell[3], stride[3], ngh[3];
  long ind[3], ip1, ip2;
  real *metric[3][2];
  real *invmet[3][2];
  real vel, slope_vel, source_vel;
  fw = CurrentFluidPatch;
  getgridsize (fw->desc, gncell, stride);
  Size = gncell[0]*gncell[1]*gncell[2];
  rho = fw->Density;
  rhop = fw->Density_Pred;
  e = fw->Energy;
  ep = fw->Energy_Pred;
  source_div = fw->SourceDiv;
  for (i = 0; i < NDIM; i++) {
    metric[i][0] = fw->desc->Metric[i][0];
    metric[i][1] = fw->desc->Metric[i][1];
    invmet[i][0] = fw->desc->InvMetric[i][0];
    invmet[i][1] = fw->desc->InvMetric[i][1];
    v[i] = fw->Velocity[i];
    vp[i] = fw->Velocity_Pred[i];
    source_v[i] = fw->SourceVelocity[i];
    slope_rho[i] = fw->SlopeDensity[i];
    slope_e[i] = fw->SlopeEnergy[i];
    for (j = 0; j < NDIM; j++) {
      slope_v[i][j] = fw->SlopeVelocity[i][j];
    }
  }
  for (i = 0; i < 3; i++)
    ngh[i] = (Nghost[i] > 0 ? Nghost[i]-1 : 0);
  for (l = 0; l < NDIM; l++) {
    for (k = ngh[2]; k < gncell[2]-ngh[2]; k++) {
      ind[2] = k;
      for (j = ngh[1]; j < gncell[1]-ngh[1]; j++) {
	ind[1] = j;
	for (i = ngh[0]; i < gncell[0]-ngh[0]; i++) {
	  ind[0] = i;
	  m = i*stride[0]+j*stride[1]+k*stride[2];
	  /* Hereafter we predict 'rho' at the cell center, at half timestep */
	  var = 0.0;
	  for (n=0; n < NDIM; n++) {
	    ip1 = (n == 0);
	    ip2 = 2 - (n == 2);
	    vel = v[n][m]; // linear or angular.
	    vel *= metric[n][0][ind[ip1]]*metric[n][1][ind[ip2]]; // converted to linear
	    slope_vel = slope_v[n][n][m]; //linear (v)-linear (dx)
	    var -= vel    * slope_rho[n][m];
	    var -= rho[m] * slope_vel;
	  }
	  var -= source_div[m]*rho[m];  //Source sign OK: evaluated as LHS.
	  rhop[m] = rho[m] + .5*dt*var;

	  /* Hereafter we predict 'v' at the cell center, at half
	     timestep, on all components */
	  for (p=0; p < NDIM; p++) { /* Loop on velocity component */
	    var = 0.0;
	    for (n=0; n < NDIM; n++) {/* Loop on slope direction */
	      ip1 = (n == 0);
	      ip2 = 2 - (n == 2);
	      vel = v[n][m]; // linear or angular.
	      vel *= metric[n][0][ind[ip1]]*metric[n][1][ind[ip2]];// converted to linear
	      slope_vel = slope_v[n][p][m];// linear-linear
	      var += -vel * slope_vel;
	    }
	    source_vel = source_v[p][m];
	    var += source_vel; //Source sign OK: evaluated as RHS.
	    /* Pressure gradient is included in source term */
	    ip1 = (p == 0);
	    ip2 = 2 - (p == 2);
	    // vp & v linear or angular. Correcting increment to match that.
	    // Will be made linear later, when filling the beam.
	    vp[p][m] = v[p][m] +.5*dt*var*invmet[p][0][ind[ip1]]*invmet[p][1][ind[ip2]];
	  }

	  /* Hereafter we predict 'internal energy' at the cell center, at half timestep */
	  var = 0.0;
	  for (n=0; n < NDIM; n++) {
	    ip1 = (n == 0);
	    ip2 = 2 - (n == 2);
	    vel = v[n][m]; // linear or angular.
	    vel *= metric[n][0][ind[ip1]]*metric[n][1][ind[ip2]]; // converted to linear
	    slope_vel = slope_v[n][n][m]; //linear (v)-linear (dx)
	    var -= vel        * slope_e[n][m];
	    var -= GetGamma()*e[m] * slope_vel;
	  }
	  var -= GetGamma()*source_div[m]*e[m];  //Source sign OK: evaluated as LHS.
	  ep[m] = e[m] + .5*dt*var;

#ifdef TEMPFLOOR
	  temperature = ep[m]/(CV*rhop[m]);
	  if (temperature < BOUNDTEMP/TEMP0/TEMPRATIO)
	    ep[m] = BOUNDTEMP/TEMP0/TEMPRATIO * CV * rhop[m];
#endif

	}
      }
    }
  }
}
