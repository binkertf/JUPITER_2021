#include "jupiter.h"

/*
 definition of the initial conditions
 (called by InitialCondition in Init.c)
*/

inline real keplerian_init(component, radius, colatitude, sigma0, a, h, f)
     int component;
     real radius, colatitude;
     real sigma0, a, h, f; 	/* sigma0, sigmaslope, aspect ratio and flaring index */
{
  real init=0.;
  real h2,  hm2, b, xi, w, isc=1., isc2f=1.;
  xi = a+1.+f;
  b = .5-f;
  /* w scales with the pressure gradient */
  w = 2.*b+xi;
  if (NDIM == 2) w= 2.*b+a;
  h2 = h*h*pow(radius,2.*f);
  hm2 = 1./h2;
  if (NDIM == 3) {
    isc = 1./sin(colatitude);
    isc2f = pow (isc, 2.*f);
  }
  switch (component) {
  case _density_:
    if (NDIM ==2)
      init = sigma0*pow(radius,-a);
    else if (NDIM ==3) {
      init = sigma0/sqrt(2.0*M_PI)/h	\
	*pow(radius,-xi) * pow(isc, w);
      if (fabs(f) < 1e-10) /* flat case */
	init *= pow(isc , -hm2);
      else		   /* flaring case */
	init *= exp(hm2 * (1. - isc2f)/2./f);
    }
    //if (fabs(init) < 1.e-18)
    //  init = 1.e-18;
    break;
  case _vazimuth_:		/* Same in 2 and 3D */
    init = pow(radius, -1.5);
    init *= sqrt(isc2f - h2*w);
    init *= isc;
    break;
  case _energy_:			/* Same in 2 and 3D */
    init = h*h*pow(radius, -2.*b);
    break;
  case _vrad_:
    break;
  	/* Viscous drift */
    init = 3.0*VISCOSITY/radius*(xi-1.5);
    if (NDIM == 2)
      init = 3.0*VISCOSITY/radius*(a-.5);
    break;
  case _vcolatitude_:
    break;
  default:
    prs_error ("Unknown component in %s at line %d", __FILE__, __LINE__);
  }
  if (component==_vazimuth_) init -= OMEGAFRAME;
  return init;
}



inline real keplerian_dust_init(component, radius, colatitude, sigma0, a, h, f)
     int component;
     real radius, colatitude;
     real sigma0, a, h, f; 	/* sigma0, sigmaslope, aspect ratio and flaring index */
{
  real init=0.;
  real h2,  hm2, b, xi, w, isc=1., isc2f=1., St_mid, alpha,hg, hd,omegakep;
  hg = h*radius;
  if (constSt==TRUE){
      St_mid = STOKESNUMBER; //3.14159265/2.0*DUSTSIZE/(sigma0/DUSTTOGAS*pow(radius,-a)); // midplane Stokes number
  }else{// constant particle size
    St_mid = M_PI/2.0*DUSTSIZE/(sigma0/DUSTTOGAS*pow(radius,-a)); // midplane Stokes number
  }

  omegakep = OMEGAFRAME*sin(colatitude)/(sqrt(radius)*sqrt(radius)*sqrt(radius));
  alpha = VISCOSITY/(omegakep*hg*hg);//*(pow(radius,0.2))*1.8;
  hd = hg*sqrt(alpha/(alpha+St_mid));
  //printf("%lg \n",hd);
  xi = a+1.+f;
  b = .5-f;
  /* w scales with the pressure gradient */
  w = 2.*b+xi;
  if (NDIM == 2) w= 2.*b+a;
  h2 = hd*hd*pow(radius,2.*f);
  hm2 = 1./h2;
  if (NDIM == 3) {
    isc = 1./sin(colatitude);
    isc2f = pow (isc, 2.*f);
  }
  switch (component) {
  case _density_:
    if (NDIM ==2)
      init = sigma0*pow(radius,-a);
    else if (NDIM ==3) {
      if(constSt==TRUE){
        init = sigma0*pow(radius,-a)/sqrt(2.0*M_PI)/hd*exp(-radius*cos(colatitude)*radius*cos(colatitude)/(2*hd*hd));
      }else{
        init = sigma0*pow(radius,-a)/sqrt(2.0*M_PI)/hd*exp(-St_mid/alpha*(exp(radius*cos(colatitude)*radius*cos(colatitude)/(2.*hg*hg))-1.0)-radius*cos(colatitude)*radius*cos(colatitude)/(2.*hg*hg)); //acc. to Fromang & Nelson 2009 Eq. (19)
      }


      //init = sigma0/sqrt(2.0*M_PI)/hd*pow(radius,-xi) * pow(isc, w);


  //    if (fabs(f) < 1e-10) /* flat case */
	//init *= pow(isc , -hm2);
  //    else		   /* flaring case */
	//init *= exp(hm2 * (1. - isc2f)/2./f);
    }


    if(init < DUSTDENSFLOOR)
      init = DUSTDENSFLOOR;

    break;
  case _vazimuth_:		/* Same in 2 and 3D */
    h2 = h*h*pow(radius,2.*f);
    init = pow(radius, -1.5);
    init *= sqrt(isc2f - h2*w);
    init *= isc;
    break;
  case _energy_:			/* Same in 2 and 3D */
    init = 0.0;
    break;
  case _vrad_:			/* Viscous drift */
    init = 3.0*VISCOSITY/radius*(xi-1.5);
    if (NDIM == 2)
      init = 3.0*VISCOSITY/radius*(a-.5);
    break;
  case _vcolatitude_:
    break;
  default:
    prs_error ("Unknown component in %s at line %d", __FILE__, __LINE__);
  }
  if (component==_vazimuth_) init -= OMEGAFRAME;
  return init;
}






/*
  gravitational external potential
  (for ComputeExternalPotential in ExtPot.c)
*/
inline real keplerian_pot(radius, azimuth, thirdcoord)
     real radius, azimuth, thirdcoord;
{
  real colatitude, z, x;
  real potential=0., planetmass, MassTaper;
  if (MASSTAPER > 1e-10) {
    MassTaper = GlobalDate/(MASSTAPER*2.0*M_PI);
    MassTaper = (MassTaper > 1.0 ? 1.0 : pow(sin(MassTaper*M_PI/2.0),2.0));
  } else
    MassTaper = 1.0;
  planetmass = PLANETMASS * MassTaper;
  switch (CoordType) {
  case CYLINDRICAL :
    z = thirdcoord;
    potential = -1./sqrt(radius*radius + z*z);
    x = radius*cos(azimuth);
    potential += -planetmass						\
      / (sqrt(1. + radius*radius -2.*radius*cos(azimuth) + z*z + SMOOTHING*SMOOTHING));
    potential += planetmass*x;
    break;
  case SPHERICAL :
    colatitude = thirdcoord;
    potential = -1./radius;
    x = radius*cos(azimuth)*sin(colatitude);
    potential += -planetmass						\
      / sqrt(1. + radius*radius -2.*radius*sin(colatitude)*cos(azimuth) + SMOOTHING*SMOOTHING);
	potential += planetmass*x;
//    / sqrt(1.000041659344951928*1.000041659344951928 + radius*radius -2.*1.000041659344951928*radius*sin(colatitude)*cos(azimuth) + SMOOTHING*SMOOTHING);
//    potential += planetmass*x*1.000041659344951928;
//printf("potential=%lg\n",potential);
    break;
  }
  return potential;
}

/*
  for a better fill of the ghosts
  (by ExecCommUp in comm_fill.c and ExecComm in comm_exec.c)
*/
inline real keplerian_comm(value, component, radius, colatitude, direction)
     real 	value;
     real 	radius, colatitude;
     int 	direction;
     long       component;
{
  real new_value = value, sc, f, h2, w, xi, d;
  d = (real)direction;
  f = FLARINGINDEX;
  if (NDIM ==3) {
    sc = sin(colatitude);
    w = 2.+SIGMASLOPE-f;	/* pressure gradient */
    xi = 1.+SIGMASLOPE+f;	/* density gradient */
  }
  else {
    sc = 1.0;
    w = 1.+SIGMASLOPE-2.*f;
    xi = SIGMASLOPE;		/* density gradient */
  }
  h2 = ASPECTRATIO*ASPECTRATIO*pow(radius, 2.*f);
  switch (component) {
  case _vazimuth_:		/* Works in 2 and 3D */
    if (direction == +1) {
      new_value += OMEGAFRAME;
      new_value *= pow(radius, 1.5)*sc;
      new_value /= sqrt(pow(sc,-2.*f)-h2*w);
    } else {
      new_value *= sqrt(pow(sc,-2.*f)-h2*w);
      new_value *= pow(radius, -1.5)/sc;
      new_value -= OMEGAFRAME;
    }
    break;
  case _energy_:		/* Works in 2 and 3D */
    new_value *= pow(radius, (real)direction*(1.-2.*f));
    break;
  case _density_:		/* Works in 2 and 3D */
    new_value *= pow(radius,d*xi)*pow(sc,d*w);
    if (fabs(f) < 1e-10)
      new_value *= pow(sc,-d/h2);
    else
      new_value *= exp(-d*(1./h2*(1.-pow(sc,-2.*f))/2./f));
    break;
  }
  return new_value;
}

/*
  better boundary conditions
  (for boundary in boundary.c)
*/
real keplerian_boundary(condition, component, value, x, xg, yg, zg)
     int condition, component;
     real value, x, xg, yg, zg;
{
  real ghost, xi, radius, cc, ccg, h2, rr;
  ghost = value;
  if (component == _density_) {
    switch (condition) {
    case 11 : //Radial boundary condition, hence : xg=ghost radius, x=active radius
      xi = SIGMASLOPE;
      if (NDIM == 3)
	xi += 1+FLARINGINDEX;
      ghost *= pow (x/xg, xi);
      break;
    case 12 : //Colatitude boundary condition
      radius = zg;
      cc = cos(x);              /* cos of colatitude of active zone */
      ccg = cos(xg);            /* cos of colatitude of ghost zone */
      h2 = ASPECTRATIO*ASPECTRATIO*pow(radius, 2.*FLARINGINDEX);
      ghost *= exp(.5*(cc*cc-ccg*ccg)/h2);	/* Approximately valid assuming a gaussian profile */
      break;
    }
  }
  if (component == _vazimuth_) {
    switch (condition) {	/* The two following prescriptions are
				   valid for a non-flaring disk. */
    case 11:
      value += OMEGAFRAME;
      ghost = value * pow(x/xg,1.5);
      ghost -= OMEGAFRAME;
      break;
    case 12:
      ghost += OMEGAFRAME;
      ghost *= sin(x)/sin(xg);
      ghost -= OMEGAFRAME;
      break;
    }
  }
  if (component == _energy_) {
    switch (condition) {
    case 11:
      ghost = value*pow(xg/x,-1.+2.*FLARINGINDEX);
      break;
    case 12:
      ghost = value;
      break;
    }
  }
  if (component == _vrad_) {
    switch (condition) {
    case 11:
      //if (x < xg) rr = x/xg;
      //else rr= xg/x;
      rr = x/xg;
      ghost = value*rr;
      break;
    case 12:
      ghost = value;
      break;
    }
  }
  if (component == _vcolatitude_) {
    switch (condition) {
    case 11:
      break;
    case 12:
      ghost *= sin(x)/sin(xg);
      break;
    }
  }
  return ghost;
}
