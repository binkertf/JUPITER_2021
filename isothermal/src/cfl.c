#include "jupiter.h"

real MaxLowLevDT[REFINEUPPERLIMIT];
boolean LevelHasChangedSinceCFL[REFINEUPPERLIMIT];
static real MaxLowLevDT_loc[REFINEUPPERLIMIT];
extern long TimeStepRatio[REFINEUPPERLIMIT];
extern long BaseStepRatio[REFINEUPPERLIMIT];


real CourantLimit (fp)
     FluidPatch *fp;
{
  long i[3], gncell[3], dim, m, b, stride[3],m_mon_crit;
  long i_mon[] = {0, 0, 0};	  /* The following variables are used only if */
  real dxmon[] = {0.0, 0.0, 0.0}; /* the -c switch is on (CFLMonitor) */
  real uxmon[] = {0.0, 0.0, 0.0};
  real xmon[] = {0.0, 0.0, 0.0};
  real cs_mon  = 0.0;
  real dxmon_crit[] = {0.0, 0.0, 0.0};
  real xmon_crit[] = {0.0, 0.0, 0.0};
  real uxmon_crit[] = {0.0, 0.0, 0.0};
  real cs_mon_crit  = 0.0;
  real dxmon_visc[] = {0.0, 0.0, 0.0};
  real xmon_visc[] = {0.0, 0.0, 0.0};
  real uxmon_visc[] = {0.0, 0.0, 0.0};
  real cs_mon_visc  = 0.0;
  long i_mon_crit[] = {0, 0, 0};
  long i_mon_visc[] = {0, 0, 0};
  boolean ViscosityLimited = NO;
  real *vel[3], *cs2, cs, dt_min=1e20, dt, sum;
  real dt_visc_min = 1e20, dt_visc = 1e20;
  real *edges[3], radius=0.0;
  real dx, u, ub, dens;
FILE *log;
FILE *log2;
//	CFLMonitor=YES;
//	printf("IN cfl %s\n",CFLMonitor);
//CFL
//
//  log2 = prs_opena ("u.dat");
  getgridsize (fp->desc, gncell, stride);
  for (dim = 0; dim < 3; dim ++) {
    vel[dim] = fp->Velocity->Field[dim];
    edges[dim] = fp->desc->Edges[dim];
  }
  cs2 = fp->Energy->Field;
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
	    //radius = fp->desc->Edges[_RAD_][i[_RAD_]];
	    radius = fp->desc->Center[_RAD_][m];
	    dx *= radius;
	    u  *= radius;
	  }
	  if (__SPHERICAL && (dim == _COLAT_)) {
	    //radius = fp->desc->Edges[_RAD_][i[_RAD_]];
	    radius = fp->desc->Center[_RAD_][m];
	    dx *= radius;
	    u  *= radius;
	  }
	  if (__SPHERICAL && (dim == _AZIM_)) {
	    //radius = fp->desc->Edges[_RAD_][i[_RAD_]];
	    radius = fp->desc->Center[_RAD_][m];
	    radius *= sin(fp->desc->Center[_COLAT_][m]);
	    dx *= radius;
	    u  *= radius;
	  }
	  if (__SPHERICAL ) {
	    radius = fp->desc->Center[_RAD_][m];
		xmon[0]=fp->desc->Center[_AZIM_][m];
		xmon[1]=radius;
		xmon[2]=fp->desc->Center[_COLAT_][m];
	  }
//	  if(CPU_Rank==3 && i[2]==18 && i[1]==2 && i[0]==14 && dim==1)
//	    fprintf(log2,"%.20g %lg %ld\n",GlobalDate,u,fp->desc->level);
	  if (Isothermal)
	    cs = sqrt(cs2[m]);
	  else
	    cs = sqrt(cs2[m]*GAMMA*(GAMMA-1.0));

    if (fp->next!=NULL)
      cs = 0.0;

	  dxmon[dim] = dx;
	  uxmon[dim] = u;
	  cs_mon = cs;
	  i_mon[dim] = i[dim];
	  sum += (fabs(u)+cs)/dx;

	  dt_visc = dx*dx/4.0/VISCOSITY;
	  if (dt_visc < dt_visc_min) {
	    dt_visc_min = dt_visc;
	    for (b = 0; b < NDIM; b++) {
	      i_mon_visc[b] = i_mon[b];
	      dxmon_visc[b] = dxmon[b];
	      uxmon_visc[b] = uxmon[b];
	      xmon_visc[b]=xmon[b];
	      cs_mon_visc   = cs_mon;
	    }
	  }
	}
	dt = 1.0/sum;
	if ((dt < dt_min) || (dt_visc_min < dt_min)) {
	  dt_min = (dt < dt_visc_min ? dt : dt_visc_min);
//	  if (CFLMonitor == YES) {
	    for (b = 0; b < NDIM; b++) {
	      dxmon_crit[b] = dxmon[b];
	      uxmon_crit[b] = uxmon[b];
	      cs_mon_crit   = cs_mon;
	      i_mon_crit[b] = i_mon[b];
              m_mon_crit=m;
	      xmon_crit[b] = xmon[b];
	    }
	    ViscosityLimited = NO;
	    if (dt_visc_min < dt) {
	      ViscosityLimited = YES;
	      for (b = 0; b < NDIM; b++) {
		dxmon_crit[b] = dxmon_visc[b];
		uxmon_crit[b] = uxmon_visc[b];
		xmon_crit[b] = xmon_visc[b];
		cs_mon_crit   = cs_mon_visc;
		i_mon_crit[b] = i_mon_visc[b];
	      }
	    }
//	  }
	}
      }
    }
  }
//  if (CFLMonitor == YES) {
	if(dt_min < 1.e-6){
    printf("-----------------------\n");
    printf ("Monitoring CFL-critical zone in patch #%ld ===> dtmin = %.12g\n",\
	    fp->desc->parent, dt_min);
    printf ("CPU_rank is : %i \n",CPU_Rank);
    printf ("Limiting zone has coordinates ( %ld  %ld  %ld m=%ld) = ( %.10g  %.10g  %.10g )\n",\
	    i_mon_crit[0], i_mon_crit[1], i_mon_crit[2], m_mon_crit,xmon_crit[0], xmon_crit[1], xmon_crit[2]);
    printf ("(active part begins at (%ld,%ld,%ld))\n", Nghost[0], Nghost[1], Nghost[2]);
    for (b = 0; b < NDIM; b++) {
      printf ("Linear vel along dim #%ld :  %.10g\n", b, uxmon_crit[b]);
      printf ("Linear zone size dim #%ld :  %.10g\n", b, dxmon_crit[b]);
    }
    printf ("Sound speed : %.10g\n", cs_mon_crit);
    if (ViscosityLimited == YES)
      printf ("Viscosity Limited\n");
//log = prs_open ("cfl.dat");
//i[2] = 38;
//    for (i[1] = 0; i[1] < gncell[1]-Nghost[1]; i[1]++) {
      /* only the active part of the mesh */
/*      for (i[0] = 0; i[0] < gncell[0]-Nghost[0]; i[0]++) {
	m = 0;
	for (dim = 0; dim < NDIM; dim++)
	  m += i[dim]*stride[dim];
	radius = fp->desc->Center[_RAD_][m];
	xmon[0]=fp->desc->Center[_AZIM_][m];
	xmon[1]=radius;
	xmon[2]=fp->desc->Center[_COLAT_][m];
	for (dim = 0; dim < NDIM; dim++) {
	  uxmon[dim]  = vel[dim][m];
	  if (dim == _AZIM_)uxmon[dim]  *= radius;
	}
       	dens= fp->Density->Field[m];
	fprintf (log, "%ld %ld %ld %ld %lg %lg %lg %lg %lg %lg %lg\n",   i[0], i[1],i[2],m,xmon[0], xmon[1], xmon[2], uxmon[0], uxmon[1], uxmon[2],dens);
	//printf (log, "%ld %ld %ld %lg %lg %lg %lg %lg %lg\n",   i[0], i[1],i[2],xmon[0], xmon[1], xmon[2], uxmon[0], uxmon[1], uxmon[2]);
}
}
 fclose (log);
 exit(1);
  }
 fclose (log2);*/
        }
  return dt_min * CFLSECURITY;
}

void UpdateCourantLimit (level)
     long level;
{
  tGrid_CPU *item;
  FluidPatch *Fluid;
  real dt;
  item = Grid_CPU_list;
  MaxLowLevDT_loc[level] = 1e20;
  while (item != NULL) {
    if (item->cpu == CPU_Rank) {
      if (level == item->level) {
	Fluid = item->Fluid;
	while (Fluid != NULL) {
	  dt = CourantLimit (Fluid);
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
//	printf("calling update for level %ld\n",i);
     UpdateCourantLimit(i);
      LevelHasChangedSinceCFL[i] = NO;
    }
  }
  MPI_Allreduce (MaxLowLevDT_loc, MaxLowLevDT, LevMax+1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  for (i = 0; i <= LevMax; i++) {
//	printf("in CLG %ld %lg %lg %ld %lg\n",i,dt_min,BaseStepRatio[i]*MaxLowLevDT[i],BaseStepRatio[i],MaxLowLevDT[i]);
    if (dt_min > BaseStepRatio[i]*MaxLowLevDT[i])
      dt_min = BaseStepRatio[i]*MaxLowLevDT[i];
//	printf("dt in CLG %lg\n",dt_min);
  }
  return (dt_min < DT ? dt_min : DT);
}
