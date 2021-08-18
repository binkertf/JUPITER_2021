#include "jupiter.h"

extern real BB_Xmin[3], BB_Xmax[3];

inline void convert_coord (x, c)
     real *x, *c;
{
  if (__CARTESIAN) {
    c[_X_] = x[0];
    c[_Y_] = x[1];
    c[_Z_] = x[2];
  }
  if (__CYLINDRICAL) {
    c[_RAD_] = sqrt(x[0]*x[0]+x[1]*x[1]);
    c[_AZIM_] = atan2(x[1], x[0]);
    c[_Vertical_] = x[2];
  }
  if (__SPHERICAL) {
    c[_RAD_] = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
    c[_AZIM_] = atan2(x[1],x[0]);
    if (c[_AZIM_] < corner_min0[_AZIM_])
      c[_AZIM_] += 2.0*M_PI;
    if (c[_AZIM_] > corner_max0[_AZIM_])
      c[_AZIM_] -= 2.0*M_PI;
    c[_COLAT_] = acos(x[2]/c[_RAD_]);
  }
}

void Ray_Integrate (Xo, Yo, Zo, ddx, ddy, ddz, ray, depth, res)
     real Xo, Yo, Zo, ddx, ddy, ddz, res;
     float *ray;
     long depth;
{
  real x[3],c[3],dx[3],i[3],s=0.0, ndx, ds;
  real xw[3], weight, wk,tau[20],kappa[20],val[200];
  real max, s_max, sum, frac_zone, dp;
  real lambda, value, *dens, tau_before_obj=0.0;
  long ii[3],ie[3],j,m,ms,mh,k,u[3],nchan,ch;
  boolean out, slice, object_found=NO;
  boolean capture_slice_eq, capture_slice_med;
  real dist_star_eq_plane, s_eq, s_med;
  real dist_planet_eq_plane, xx, yy, zz;
  real slice_med=0.0, slice_eq=0.0;
  tGrid_CPU *item;
  tGrid *g;
  frac_zone = .1;
  if (CoarseRayTracing) frac_zone = .8;
  x[0] = Xo;
  x[1] = Yo;
  x[2] = Zo;
  ndx = sqrt(ddx*ddx+ddy*ddy+ddz*ddz);
  ds = (RANGE1HIGH-RANGE1LOW)/(real)SIZE1*0.5;
  ddx /= ndx;
  ddy /= ndx;
  ddz /= ndx;
  if (ddz == 0.0) {
    dist_star_eq_plane = -1.0;
    dist_planet_eq_plane = -1.0;
    s_eq = -1.0;
  } else {
    xx = Xo+ddx*(0.-Zo)/ddz;
    yy = Yo+ddy*(0.-Zo)/ddz;
    dist_planet_eq_plane = sqrt((xx-1.0)*(xx-1.0)+yy*yy);
    dist_star_eq_plane = sqrt(xx*xx+yy*yy);
    s_eq = sqrt((xx-Xo)*(xx-Xo)+(yy-Yo)*(yy-Yo)+Zo*Zo);
    if (ddz*Zo > 0.0) s_eq = -s_eq;
  }
  if (ddy == 0.0) {
    s_med = -1.0;
  } else {
    xx = Xo+ddx*(0.-Yo)/ddy;
    zz = Zo+ddz*(0.-Zo)/ddz;
    s_med = sqrt((xx-Xo)*(xx-Xo)+Yo*Yo+(zz-Zo)*(zz-Zo));
    if (ddy*Yo > 0.0) s_med = -s_med;
  }
  dx[0] = ds*ddx;
  dx[1] = ds*ddy;
  dx[2] = ds*ddz;
  make_background_stars (dx, ray, res);
  max = 0.0;
  sum = 0.0;
  s_max = 250.0;
  nchan = depth-9;
  if (nchan > 200) prs_error ("Too many opacity channels\n");
  for (ch = 0; ch < nchan; ch++) {
    kappa[ch] = 1./SIGMA0*(ASPECTRATIO+2.*(real)ch/(real)nchan);
    val[ch] = tau[ch] = 0.0;
  }
  while (s < s_max) {
    out = FALSE;
    xw[0] = x[0];
    xw[1] = x[1];
    xw[2] = x[2];
    if (RayMirror) {
      if (x[2] < 0.0) xw[2] = -x[2];
    }
    for (j = 0; j < 3; j++) {
      if (((x[j] < BB_Xmin[j]) && dx[j] < 0.0) ||	\
	  ((x[j] > BB_Xmax[j]) && dx[j] > 0.0)) {
	out = TRUE;
	s = s_max;
      }
    }
    for (j = 0; j < 3; j++) {
      if ((x[j] < BB_Xmin[j]) && (dx[j] > 0.0)) {
	lambda = (BB_Xmin[j]-x[j])/dx[j]+1e-9;
	for (k = 0; k < 3; k++)
	  x[k] += lambda*dx[k];
	s += ds*lambda;
      } else if ((x[j] > BB_Xmax[j]) && (dx[j] < 0.0)) {
	lambda = (BB_Xmax[j]-x[j])/dx[j]+1e-9;
	for (k = 0; k < 3; k++)
	  x[k] += lambda*dx[k];
	s += ds*lambda;
      }
    }
    if (!out) {
      slice = NO;
      convert_coord (xw,c);
      dp = (sqrt((x[0]-1.0)*(x[0]-1.0)+x[1]*x[1]+x[2]*x[2])*1.14);
      if (dp < 2e-4) {
	if (!object_found)
	  make_planet (dx, x, c, ray);
	object_found = YES;
      }
      if (c[_RAD_] < .1) {
	if (!object_found)
	  make_central_star (dx, x, c, ray);
	object_found = YES;
      }
      item = Grid_CPU_list;
      while (item != NULL) {
	/* Are we in the CPU grid 'item' ? */
	/* Note : we scan the CPU grids instead of the grids, */
	/* despite the fact that we know that CPU_Number = 1, */
	/* because the 'overlap' flag is set for CPU_grids */
	g = item->Parent;
	if ((c[0] >= g->corner_min[0]) && (c[0] <= g->corner_max[0]) &&	\
	    (c[1] >= g->corner_min[1]) && (c[1] <= g->corner_max[1]) &&	\
	    (c[2] >= g->corner_min[2]) && (c[2] <= g->corner_max[2])) {
	  dens = item->Fluid->Density->Field;
	  for (j = 0; j < 3; j++) {
	    i[j] = (real)NGH+(c[j]-g->corner_min[j])/\
	      (g->corner_max[j]-g->corner_min[j]+1e-7)*g->ncell[j];
	    ie[j] = (long)(i[j]);
	    i[j] -= .5;
	    ii[j] = (long)(i[j]);
	  }
	  m = mh = 0;
	  for (j = 0; j < 3; j++) {
	    mh += ie[j]*item->stride[j];
	    m += ii[j]*item->stride[j];
	  }
	  if (!(item->Hidden[mh])) {
	    ds = frac_zone*(item->Edges[_RAD_][2]-item->Edges[_RAD_][1]);
	    capture_slice_med = NO;
	    if ((s < s_med) && (s+ds >= s_med)) {
	      ds=s_med-s+1e-5;
	      capture_slice_med = YES;
	    }
	    capture_slice_eq = NO;
	    if ((s < s_eq) && (s+ds >= s_eq)) {
	      ds=s_eq-s+1e-5;
	      capture_slice_eq = YES;
	      capture_slice_med = NO;
	    }
	  /* The above expression could be refined */
	    dx[0] = ds*ddx;
	    dx[1] = ds*ddy;
	    dx[2] = ds*ddz;
	    value = 0.0;
	    for (u[2] = 0; u[2] < 2; u[2]++) {
	      for (u[1] = 0; u[1] < 2; u[1]++) {
		for (u[0] = 0; u[0] < 2; u[0]++) {
		  weight = 1.0;
		  ms = m;
		  for (k = 0; k < 3; k++) {
		    wk = i[k]-(real)ii[k];
		    weight *= (u[k] == 0 ? 1.0-wk : wk);
		    if (u[k] > 0) ms += item->stride[k];
		  }
		  value += weight*dens[ms];
		}
	      }
	    }
	    if (capture_slice_med) slice_med = value;
	    if (capture_slice_eq) slice_eq = value;
	    sum += value*ds;
	    for (ch = 0; ch < nchan; ch++) {
	      val[ch] += ds*value*exp(-tau[ch]);
	      tau[ch] += ds*value*kappa[ch];
	    }
	    if (!object_found) tau_before_obj = tau[0];
	    if (value > max)
	      max = value;
	  }
	}
	item = item->next;
      }
    }
    for (j = 0; j < 3; j++)
      x[j] += dx[j];
    s += ds;
  }
  ray[2] = (float)max;
  ray[3] = (float)sum;
  ray[4] = tau_before_obj;
  ray[5] = slice_med;
  ray[6] = slice_eq;
  ray[7] = dist_star_eq_plane;
  ray[8] = dist_planet_eq_plane;
  for (ch = 0; ch < nchan; ch++) {
    ray[9+ch] = (float)val[ch];
  }
}
