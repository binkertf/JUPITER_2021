#include "jupiter.h"

#define NBSTARS 6000

static real st_col[NBSTARS], st_phi[NBSTARS], st_lum[NBSTARS];

void make_stars ()
{
  long i;
  real mag, ct;
  srand48 (56);
  for (i = 0; i < NBSTARS; i++) {
    ct = drand48()*2.-1.;
    st_col[i] = acos (ct);
    st_phi[i] = drand48()*2.*M_PI;
    mag = -2.0;
    while (mag < -1.0) {
      mag = 6.0+.9*log(drand48());
    }
    st_lum[i] = pow(2.512,-mag);
  }
}

void make_central_star (dx, x, c, ray)
     real *dx, *x, *c;
     float *ray;
{
  double a1=.93, a2=-.23;
  double xc[3] = {0., 0., 0.}, xn[3];
  double scal;
  (void) c;
  xn[0] = x[0] - xc[0];
  xn[1] = x[1] - xc[1];
  xn[2] = x[2] - xc[2];
  scal = -dx[0]*xn[0]-dx[1]*xn[1]-dx[2]*xn[2];
  scal /= sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  scal /= sqrt(xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]);
  ray[1] = (float)(1.-a1-a2+a1*scal+a2*scal*scal);
}

void make_planet (dx, x, c, ray)
     real *dx, *x, *c;
     float *ray;
{
  double a1=.93, a2=-.23;
  double xc[3] = {1., 0., 0.}, xn[3];
  double scal;
  (void) c;
  xn[0] = x[0] - xc[0];
  xn[1] = x[1] - xc[1];
  xn[2] = x[2] - xc[2];
  scal = -dx[0]*xn[0]-dx[1]*xn[1]-dx[2]*xn[2];
  scal /= sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  scal /= sqrt(xn[0]*xn[0]+xn[1]*xn[1]+xn[2]*xn[2]);
  ray[1] = (float)(1.-a1-a2+a1*scal+a2*scal*scal);
}

void make_background_stars (dx, ray, res)
     real *dx, res;
     float *ray;
{
  real col, phi, r, d_phi, d_col;
  static boolean FirstTime = YES;
  long i;
  ray[0] = 0.0;
  if (FirstTime) make_stars ();
  FirstTime = NO;
  r = sqrt(dx[0]*dx[0]+dx[1]*dx[1]+dx[2]*dx[2]);
  phi = atan2(dx[1], dx[0]);
  if (phi < 0.0) phi+=2.0*M_PI;
  col = acos(dx[2]/r);
  for (i = 0; i < NBSTARS; i++) {
    d_phi = phi-st_phi[i];
    d_col = col-st_col[i];
    if ((fabs(d_phi) < 2.*res) &&		\
	(fabs(d_col) < 2.*res)) {
      d_phi *= sin(col);
      r = sqrt(d_phi*d_phi+d_col*d_col);
      ray[0] += st_lum[i]*exp(-r*r/res/res*2.0);
    }
  }
}
