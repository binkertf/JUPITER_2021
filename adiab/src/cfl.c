#include "jupiter.h"

real MaxLowLevDT[REFINEUPPERLIMIT];
boolean LevelHasChangedSinceCFL[REFINEUPPERLIMIT];
static real MaxLowLevDT_loc[REFINEUPPERLIMIT];
extern long TimeStepRatio[REFINEUPPERLIMIT];
extern long BaseStepRatio[REFINEUPPERLIMIT];


real CourantLimit (fp)
     FluidPatch *fp;
{
  long i[3], gncell[3], dim, m, b, stride[3];
  long i_mon[] = {0, 0, 0};	  /* The following variables are used only if */
  real dxmon[] = {0.0, 0.0, 0.0}; /* the -c switch is on (CFLMonitor) */
  real uxmon[] = {0.0, 0.0, 0.0};
  real cs_mon=0.0, rho_mon=0.0, ene_mon=0.0;
  real dxmon_crit[] = {0.0, 0.0, 0.0};
  real uxmon_crit[] = {0.0, 0.0, 0.0};
  real cs_mon_crit  = 0.0, rho_mon_crit=0.0, ene_mon_crit=0.0;
  real dxmon_visc[] = {0.0, 0.0, 0.0};
  real uxmon_visc[] = {0.0, 0.0, 0.0};
  real cs_mon_visc  = 0.0;
  long i_mon_crit[] = {0, 0, 0};
  long i_mon_visc[] = {0, 0, 0};
  boolean ViscosityLimited = NO;
  real *vel[3], *cs2, cs, *rho, *temp, *gamma, dt_min=1e20, dt, sum;
  real dt_visc_min = 1e20, dt_visc = 1e20;
  real *edges[3], radius=0.0;
  real dx, u;
  getgridsize (fp->desc, gncell, stride);
  for (dim = 0; dim < 3; dim ++) {
    vel[dim] = fp->Velocity->Field[dim];
    edges[dim] = fp->desc->Edges[dim];
  }
  cs2 = fp->Energy->Field;
  rho = fp->Density->Field;
  temp = fp->Temperature->Field;
  //gamma = fp->Gamma->Field;
  for (i[2] = Nghost[2]; i[2] < gncell[2]-Nghost[2]; i[2]++) {
    /* We do not scan the ghosts */
    for (i[1] = Nghost[1]; i[1] < gncell[1]-Nghost[1]; i[1]++) {
      /* only the active part of the mesh */
      for (i[0] = Nghost[0]; i[0] < gncell[0]-Nghost[0]; i[0]++) {
	  	m = 0;
		for (dim = 0; dim < NDIM; dim++)
	  		m += i[dim]*stride[dim];
		sum = 0.0;
		for (dim = 0; dim < NDIM; dim++) {
	  		dx = edges[dim][i[dim]+1]-edges[dim][i[dim]];
	  		u  = vel[dim][m];
	  	if (__CYLINDRICAL && (dim == _AZIM_)) {
	    	//radius = fp->desc->Edge[_RAD_][i[_RAD_]];
	    	radius = fp->desc->Center[_RAD_][m];
	    	dx *= radius;
	    	u  *= radius;
	  	}
	  	if (__SPHERICAL && (dim == _COLAT_)) {
	    	//radius = fp->desc->Edge[_RAD_][i[_RAD_]];
	    	radius = fp->desc->Center[_RAD_][m];
	    	dx *= radius;
	    	u  *= radius;
	  	}
	  	if (__SPHERICAL && (dim == _AZIM_)) {
	    	//radius = fp->desc->Edge[_RAD_][i[_RAD_]];
	    	radius = fp->desc->Center[_RAD_][m];
	    	radius *= sin(fp->desc->Center[_COLAT_][m]);
	    	dx *= radius;
	    	u  *= radius;
	  	}
	  	if (Isothermal || (strncasecmp(fp->Name, "dust", 4) == 0))
	    	cs = sqrt(cs2[m]); //isothermal sound speed
	  	else{
			// Gamma = Gamma1(double temp, double rho) -> temp[m], rho[m] are both of type real 
			cs = sqrt(cs2[m]/rho[m]*GetGamma()*(GetGamma()-1.0)); //adiabatic sound speed
	  	}
	  dxmon[dim] = dx;
	  uxmon[dim] = u;
	  cs_mon = cs;
	  rho_mon = rho[m];
	  ene_mon = cs2[m];
	  i_mon[dim] = i[dim];
	  if ((strncasecmp(fp->Name, "dust", 4) == 0) && (diffmode != 1)){  
		sum += (fabs(u))/dx;
	  }else{
		sum += (fabs(u)+cs)/dx;
	  }
	
	  dt_visc = dx*dx/4.0/VISCOSITY;
	  if (dt_visc < dt_visc_min) {
	    dt_visc_min = dt_visc;
	    for (b = 0; b < NDIM; b++) {
	      i_mon_visc[b] = i_mon[b];
	      dxmon_visc[b] = dxmon[b];
	      uxmon_visc[b] = uxmon[b];
	      cs_mon_visc   = cs_mon;
	    }
	  }
	}
	dt = 1.0/sum;
	if ((dt < dt_min) || (dt_visc_min < dt_min)) {
	  dt_min = (dt < dt_visc_min ? dt : dt_visc_min);
	  if (CFLMonitor == YES) {
	    for (b = 0; b < NDIM; b++) {
	      dxmon_crit[b] = dxmon[b];
	      uxmon_crit[b] = uxmon[b];
	      cs_mon_crit   = cs_mon;
	      rho_mon_crit  = rho_mon;
	      ene_mon_crit  = ene_mon;
	      i_mon_crit[b] = i_mon[b];
	    }
	    ViscosityLimited = NO;
	    if (dt_visc_min < dt) {
	      ViscosityLimited = YES;
	      for (b = 0; b < NDIM; b++) {
		dxmon_crit[b] = dxmon_visc[b];
		uxmon_crit[b] = uxmon_visc[b];
		cs_mon_crit   = cs_mon_visc;
		i_mon_crit[b] = i_mon_visc[b];
	      }
	    }
	  }
	}
      }
    }
  }
  if (CFLMonitor == YES) {
    printf ("Monitoring CFL-critical zone in patch #%ld ===> dtmin = %.12g\n",\
	    fp->desc->parent, dt_min);
    printf ("Limiting zone has coordinates ( %ld  %ld  %ld )\n",\
	    i_mon_crit[0], i_mon_crit[1], i_mon_crit[2]);
    printf ("(active part begins at (%ld,%ld,%ld))\n", Nghost[0], Nghost[1], Nghost[2]);
    for (b = 0; b < NDIM; b++) {
      printf ("Linear vel along dim #%ld :  %.10g\n", b, uxmon_crit[b]);
      printf ("Linear zone size dim #%ld :  %.10g\n", b, dxmon_crit[b]);
    }
    printf ("Sound speed : %.10g\n", cs_mon_crit );
    if (ViscosityLimited == YES)
      printf ("Viscosity Limited\n");
    printf ("Energy : %.10g\n", ene_mon_crit);
    printf ("Density : %.10g\n",rho_mon_crit);
  }
  return dt_min * CFLSECURITY;
}

real StoppingTimeLimit (fp)
     FluidPatch *fp;
{
  long i[3], gncell[3], dim, m, b, stride[3];
  boolean ViscosityLimited = NO;
  real *vel[3], *cs2, cs, *rho, tau_s_min=1e20, dt, sum, tau_s, dustsz, dustsolidrho, omegakep;
  real *edges[3], radius=0.0;
  getgridsize (fp->desc, gncell, stride);

  dustsz = DUSTSIZE / R0; // diameter of dust grains in code units
  dustsolidrho = DUSTSOLIDRHO / RHO0; // solid density of dust grains in code units, typically 3 g/cm^3 in physical units


  for (dim = 0; dim < 3; dim ++) {
    vel[dim] = fp->Velocity->Field[dim];
    edges[dim] = fp->desc->Edges[dim];
  }
  cs2 = fp->Energy->Field;
  rho = fp->Density->Field;
  for (i[2] = Nghost[2]; i[2] < gncell[2]-Nghost[2]; i[2]++) {
    /* We do not scan the ghosts */
    for (i[1] = Nghost[1]; i[1] < gncell[1]-Nghost[1]; i[1]++) {
      /* only the active part of the mesh */
      for (i[0] = Nghost[0]; i[0] < gncell[0]-Nghost[0]; i[0]++) {
	  	m = 0;
		for (dim = 0; dim < NDIM; dim++)
	  		m += i[dim]*stride[dim];
		radius = fp->desc->Center[_RAD_][m];
		tau_s = 0.0;
		
			if (Isothermal){//isothermal
	    		cs = sqrt(cs2[m]); //isothermal sound speed
				if (NDIM ==1){
					tau_s = STOKESNUMBER;
				}
				if (NDIM ==2){
					omegakep = 1.0/sqrt(radius*radius*radius);
					if(constSt==TRUE){
						tau_s = STOKESNUMBER / omegakep;
					}else{
						tau_s = M_PI / 2.0 * dustsz * dustsolidrho / rho[m] / omegakep;
					}
				}
				if (NDIM ==3){
					if(constSt==TRUE){//constant Stokes number
						omegakep = 1.0/sqrt(radius*radius*radius);
						tau_s = STOKESNUMBER / omegakep;
					}else{
						tau_s = sqrt(M_PI / 8.0) * dustsz * dustsolidrho / (cs * rho[m]);
					}

				}
			}else{//radiative
				// Gamma = Gamma1(double temp, double rho) -> temp[m], rho[m] are both of type real
	    		cs = sqrt(cs2[m]/rho[m]*GetGamma()*(GetGamma()-1.0)); //adiabatic sound speed
				tau_s = sqrt(GetGamma() * M_PI / 8.0) * dustsz * dustsolidrho / (cs * rho[m]);
			}
	  	if (tau_s < tau_s_min){
		  tau_s_min = tau_s;
		}

	  }
    }
  }
  return tau_s_min * CFLSECURITY;
}

void UpdateCourantLimit (level)
     long level;
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  real dt, t_stop_min;
  item = Grid_CPU_list;
  MaxLowLevDT_loc[level] = 1e20;
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      if (level == item->level) {
	Fluid = item->Fluid;
	while (Fluid != NULL) {
	  dt = CourantLimit (Fluid);
	  // gas-dust stopping time
	  if ((Fluid->next==NULL) && (SAMPLETSTOP==YES) && (NbFluids > 1)){ //gas
	  	t_stop_min = StoppingTimeLimit (Fluid); //minimum stopping time
		  
		if (t_stop_min<dt){
			dt = t_stop_min;
		}

	  }
	  if (dt < MaxLowLevDT_loc[level]) MaxLowLevDT_loc[level] = dt;
	  Fluid = Fluid->next;
	}
      }
    }
    item = item->next;
  } /* Note that at this stage the MaxLowLevDT depends on the p.e. */
}

real CourantLimitGlobal ()
{
  long i;
  real dt_min = 1e20;
  for (i = 0; i <= LevMax; i++) {
    if (LevelHasChangedSinceCFL[i] == YES) {
      UpdateCourantLimit(i);
      LevelHasChangedSinceCFL[i] = NO;
    }
  }
  MPI_Allreduce (MaxLowLevDT_loc, MaxLowLevDT, LevMax+1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  for (i = 0; i <= LevMax; i++) {
    if (dt_min > BaseStepRatio[i]*MaxLowLevDT[i])
      dt_min = BaseStepRatio[i]*MaxLowLevDT[i];
  }
  return (dt_min < DT ? dt_min : DT);
}
