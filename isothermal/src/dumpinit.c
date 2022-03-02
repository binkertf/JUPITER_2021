/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "iniparser.pl" for details   */
/*                              */
/********************************/

#include "init.h"
#include "initcodes.h"
#include "jupiter.h"

void DumpInitCode (filename)
char *filename;
{
    FILE *init;
    long i;
    init = prs_open (filename);
    for (i = 0; i < NbFluids; i++) {
    findcodeinit (InitCodeNames[i], &InitCode);
    fprintf (init, "Fluid %ld (%s):\n", i, FluidName[i]);
    switch (InitCode) {
case multi1d:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   if (x < 2.5) density = 1.0;\n");
	fprintf (init, "   if (x >= 2.5) density = 0.5;\n");
	fprintf (init, "   if (x < 2.5) vx = 0.5;\n");
	fprintf (init, "   if (x >= 2.5) vx = 0.25;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case multi1db:
	fprintf (init, "   density = 1.0e-10;\n");
	fprintf (init, "   if (x > 0.5 && x < 2.0) density = 1.0;\n");
	fprintf (init, "   if (x > 2.5 && x < 4.0) density = 0.5;\n");
	fprintf (init, "   if (x > 0.5 && x < 2.0) vx = 1.0;\n");
	fprintf (init, "   if (x > 2.5 && x < 4.0) vx = -0.2;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case multi1db_dust:
	fprintf (init, "   density = 1.0e-10;\n");
	fprintf (init, "   if (x > 0.5 && x < 2.0) density = 0.9;\n");
	fprintf (init, "   if (x > 2.5 && x < 4.0) density = 0.4;\n");
	fprintf (init, "   if (x > 0.5 && x < 2.0) vx = 1.0;\n");
	fprintf (init, "   if (x > 2.5 && x < 4.0) vx = -0.2;\n");
	fprintf (init, "   energy = 0.0;\n");
	break;
case W19_1_g:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   vx = 1.0;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case W19_1_d:
	fprintf (init, "   density = 1.0+0.001*cos(x*6.28318530718);\n");
	fprintf (init, "   vx = 1.0;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case W19_2_g:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   vx = 1.0;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case W19_2_d:
	fprintf (init, "   density = 1.0+6.451433001183836*0.001*sin(x*6.28318530718);\n");
	fprintf (init, "   vx = 1.0+0.001*cos(x*6.28318530718);\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case sod1d:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   if (x > 3.5) density = 0.1;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case sod1db:
	fprintf (init, "   density = 0.8;\n");
	fprintf (init, "   if (x > 5.5) density = 0.3;\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case aw1d:
	fprintf (init, "   density = 1.0+0.01*cos(2.0*M_PI*x/10.0);\n");
	fprintf (init, "   vx      = 0.01*cos(2.0*M_PI*x/10.0);\n");
	fprintf (init, "   energy  = CS*CS;\n");
	break;
case sphaw1d:
	fprintf (init, "   energy = CS*CS;\n");
	fprintf (init, "   density	= 1.0+0.001*exp(-(radius-5.5)*(radius-5.5)*10.0);\n");
	break;
case daw1d:
	fprintf (init, "   energy = CS*CS;\n");
	fprintf (init, "   density = 1.0+0.001*cos(2.0*M_PI*x/10.0);\n");
	break;
case awsph1dlat:
	fprintf (init, "   density = 1.0+0.001*exp(-(colatitude-M_PI*.5)*\\n");
	fprintf (init, "   (colatitude-M_PI*.5)*200.0);\n");
	fprintf (init, "   energy  = CS*CS;\n");
	break;
case vkring1d:
	fprintf (init, "   density = exp(-(radius-1.0)*(radius-1.0)/0.012)/\\n");
	fprintf (init, "   2.0/M_PI/sqrt(M_PI*0.012)/radius*sqrt(sqrt(radius));\n");
	fprintf (init, "   energy = CS*CS;\n");
	fprintf (init, "   vrad = -1.5*VISCOSITY/radius+6.0*VISCOSITY*(radius-1.0)/0.012;\n");
	fprintf (init, "   v_azimuth = 1./radius/sqrt(radius)-OMEGAFRAME;\n");
	break;
case soundwave1:
	fprintf (init, "   interm1 = cos(2.0*M_PI*x);\n");
	fprintf (init, "   density = 2.+1e-4*interm1;\n");
	fprintf (init, "   energy = 1.0;\n");
	fprintf (init, "   vx = 5e-5*interm1;\n");
	break;
case soundwave2:
	fprintf (init, "   interm1 = cos(2.0*M_PI*x);\n");
	fprintf (init, "   density = 1.-1e-4*interm1;\n");
	fprintf (init, "   energy = 1.0;\n");
	fprintf (init, "   vx = -1e-4*interm1;\n");
	break;
case sodpm06:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   vx = 0.25;\n");
	fprintf (init, "   energy = CS*CS;\n");
	fprintf (init, "   if (x > 0.0) {\n");
	fprintf (init, "   density = 0.5;\n");
	fprintf (init, "   vx = 0.5;\n");
	fprintf (init, "   }\n");
	break;
case dustpm06:
	fprintf (init, "   density   = keplerian_init(_density_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, 0.000001, 0.0);\n");
	fprintf (init, "   energy    = keplerian_init(_energy_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, 0.000001, 0.0);\n");
	fprintf (init, "   vrad      = keplerian_init(_vrad_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, 0.000001, 0.0);\n");
	fprintf (init, "   v_azimuth = keplerian_init(_vazimuth_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, 0.000001, 0.0);\n");
	fprintf (init, "   v_colatitude= keplerian_init(_vcolatitude_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, 0.000001, 0.0);\n");
	break;
case tiltawsph2d:
	fprintf (init, "   x = cos(azimuth)*sin(colatitude)*sin(0.4)+cos(colatitude)*cos(0.4);\n");
	fprintf (init, "   /* Since we are in spherical coordinates, */\n");
	fprintf (init, "   /* x can be used as a work variable */\n");
	fprintf (init, "   density = 1.0+0.001*exp(-x*x*50.0);\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case vkoffring2d:
	fprintf (init, "   energy = CS*CS;\n");
	fprintf (init, "   z = radius*cos(colatitude);\n");
	fprintf (init, "   axial_radius = radius*sin(colatitude);\n");
	fprintf (init, "   v_azimuth = pow(axial_radius, -3.0);\n");
	fprintf (init, "   v_azimuth -= energy/axial_radius*2.0*(axial_radius-1.0)/0.018;\n");
	fprintf (init, "   v_azimuth = sqrt(v_azimuth) - OMEGAFRAME;\n");
	fprintf (init, "   density = .5/M_PI/sqrt(M_PI*0.018)*sqrt(sqrt(axial_radius))/axial_radius;\n");
	fprintf (init, "   density*= exp(-pow(axial_radius-1.0,2.0)/0.018)\\n");
	fprintf (init, "   *exp(-(z-1.0)*(z-1.0)/2./.15);\n");
	fprintf (init, "   vx = -1.5*VISCOSITY/axial_radius+6.0*VISCOSITY*(axial_radius-1.0)/0.018;\n");
	fprintf (init, "   vrad = vx*sin(colatitude);\n");
	fprintf (init, "   v_colatitude = vx/radius*cos(colatitude);\n");
	break;
case sod2d:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   if (x/10.0+y/6.0 > 1.0) density = 0.1;\n");
	fprintf (init, "   density *= (1.+.0*drand48());\n");
	fprintf (init, "   energy = CS*CS;\n");
	break;
case kepler3d:
	fprintf (init, "   density   = keplerian_init(_density_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   energy    = keplerian_init(_energy_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   vrad      = keplerian_init(_vrad_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   v_azimuth = keplerian_init(_vazimuth_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   v_colatitude= keplerian_init(_vcolatitude_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	break;
case keplerdust:
	fprintf (init, "   density   = keplerian_dust_init(_density_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   energy    = 0.0;\n");
	fprintf (init, "   vrad      = 0.0;\n");
	fprintf (init, "   v_azimuth = keplerian_dust_init(_vazimuth_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   v_colatitude= 0.0; //keplerian_dust_init(_vcolatitude_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	break;
case keplerdust_2:
	fprintf (init, "   density   = keplerian_dust_init(_density_,radius,colatitude, SIGMA0*DUSTTOGAS, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   energy    = keplerian_dust_init(_energy_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   vrad      = 0.0; //keplerian_init(_vrad_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   v_azimuth = keplerian_dust_init(_vazimuth_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   v_colatitude= 0.0; //keplerian_init(_vcolatitude_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	break;
case essai:
	fprintf (init, "   density   = SIGMA0;\n");
	fprintf (init, "   vrad      = 0.0;\n");
	fprintf (init, "   energy    = ASPECTRATIO*ASPECTRATIO;\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)-OMEGAFRAME;\n");
	break;
case tides:
	fprintf (init, "   density = exp(-pow(radius,1.5));\n");
	fprintf (init, "   energy = 1./radius/radius;\n");
	break;
case homogeneity1:
	fprintf (init, "   density   = keplerian_init(_density_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   energy    = ASPECTRATIO*ASPECTRATIO/radius;\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)-OMEGAFRAME;\n");
	break;
case homogeneity2:
	fprintf (init, "   density   = keplerian_init(_density_,radius,colatitude, SIGMA0, SIGMASLOPE, ASPECTRATIO, FLARINGINDEX);\n");
	fprintf (init, "   energy    = ASPECTRATIO*ASPECTRATIO/radius*1.5e11*\\n");
	fprintf (init, "   pow((1.5e11/5.02989563261132871975e6),2.0);\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)*pow(1.5e11,1.5)/\\n");
	fprintf (init, "   5.02989563261132871975e6-OMEGAFRAME;\n");
	break;
case flatatm:
	fprintf (init, "   density   = 1.0;\n");
	fprintf (init, "   energy    = 1.0;\n");
	break;
case flatatmpl:
	fprintf (init, "   density   = 1.0+1e-4*exp(-(z-1.0)*(z-1.0)*40.0);\n");
	fprintf (init, "   energy    = 1.0;\n");
	break;
case isoatm:
	fprintf (init, "   density   = exp(-z);\n");
	fprintf (init, "   if (z > 2.0) density = exp (z-4);\n");
	fprintf (init, "   energy    = 1.0;\n");
	break;
case isoatmpert:
	fprintf (init, "   density   = exp(-z)+1e-13*exp(-(z-1.)*(z-1.)*40.);\n");
	fprintf (init, "   if (z > 2)\n");
	fprintf (init, "   density = exp(z-4)+1e-13*exp(-(z-1.)*(z-1.)*40.);\n");
	fprintf (init, "   energy    = 1.0;\n");
	break;
case isoatmpl:
	fprintf (init, "   density   = 1e-4*exp(-(z-1.0)*(z-1.0)*40.0);\n");
	fprintf (init, "   energy    = 0.0;\n");
	break;
case isoatmpll:
	fprintf (init, "   density   = 1e-2*exp(-(z-1.0)*(z-1.0)*40.0);\n");
	fprintf (init, "   energy    = 0.0;\n");
	break;
case isoatmps:
	fprintf (init, "   density   = 1e-13*exp(-(z-1.0)*(z-1.0)*40.0);\n");
	fprintf (init, "   energy    = 0.0;\n");
	break;
case isoatm2d:
	fprintf (init, "   density   = exp(-(2.*x+y));\n");
	fprintf (init, "   energy    = 1.0;\n");
	break;
case isoatm2dp:
	fprintf (init, "   density   = 1e-10*exp(-((x-1.0)*(x-1.0)+(y-1.0)*(y-1.0))*40.0);\n");
	fprintf (init, "   energy    = 0.0;\n");
	break;
case trivialsphere:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   energy = 1.0;\n");
	fprintf (init, "   v_azimuth = -OMEGAFRAME;\n");
	break;
case HVCCyl:
	fprintf (init, "   density = 2.24*exp (-z*z/2.0/70./70.);\n");
	fprintf (init, "   v_azimuth = 200.0/radius-OMEGAFRAME;\n");
	fprintf (init, "   energy = 100.0;\n");
	fprintf (init, "   x = pow((radius-7000.0),2.0)+pow(z-250.0,2.0)+pow((radius*azimuth-250.0),2.0);\n");
	fprintf (init, "   x = sqrt (x);\n");
	fprintf (init, "   if (x < 100.0) {\n");
	fprintf (init, "   density += 0.224*exp(-x*x/50./50.);\n");
	fprintf (init, "   vz   = -71.0;\n");
	fprintf (init, "   v_azimuth = v_azimuth-71.0/radius;\n");
	fprintf (init, "   }\n");
	break;
case testrot:
	fprintf (init, "   density = 1.0/radius;\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO/radius;\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)*sqrt(1.-2.*ASPECTRATIO*ASPECTRATIO)-OMEGAFRAME;\n");
	break;
case testrefleq:
	fprintf (init, "   density = exp(-(z-M_PI/2.0)*(z-M_PI/2.0)/2.0/ASPECTRATIO/ASPECTRATIO);\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO/radius;\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO;\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)*sqrt(1.-2.*ASPECTRATIO*ASPECTRATIO)-OMEGAFRAME;\n");
	fprintf (init, "   v_colatitude = 0.0*1e-10*(z-M_PI/2.0);\n");
	break;
case plan1d:
	fprintf (init, "   density = pow(z,-SMOOTHING);\n");
	fprintf (init, "   energy = 1.0;\n");
	break;
case star1d:
	fprintf (init, "   density = exp(-SMOOTHING*z*z/2.0);\n");
	fprintf (init, "   energy = 1.0;\n");
	break;
case local:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO;\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)-OMEGAFRAME;\n");
	break;
case localgrad:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO/radius;\n");
	fprintf (init, "   v_azimuth = (1.0/radius/sqrt(radius))*sqrt(1.0-ASPECTRATIO*ASPECTRATIO)-OMEGAFRAME;\n");
	break;
case novortgrad:
	fprintf (init, "   density = 1.0/radius/sqrt(radius);\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO*radius*sqrt(radius);\n");
	fprintf (init, "   v_azimuth = 1.0/radius/sqrt(radius)-OMEGAFRAME;\n");
	break;
case sod2d_vert:
	fprintf (init, "   energy = 1.0;\n");
	fprintf (init, "   density = 0.001;\n");
	fprintf (init, "   if (x < 0.5) density = 100.0;\n");
	break;
case throop:
	fprintf (init, "   radius = sqrt(x*x+y*y+z*z);\n");
	fprintf (init, "   colatitude = acos((x*INCLINATION+z)/(sqrt(1.+INCLINATION*INCLINATION)*radius));\n");
	fprintf (init, "   interm2 = cos(colatitude)*cos(colatitude)/2./ASPECTRATIO/ASPECTRATIO;\n");
	fprintf (init, "   if (interm2 > 17.0) interm2 = 17.0;\n");
	fprintf (init, "   interm1 = exp(-interm2);\n");
	fprintf (init, "   v_azimuth=1./sqrt(radius*radius*radius+1e-2)/sqrt(1.+INCLINATION*INCLINATION);\n");
	fprintf (init, "   vx = -v_azimuth*y;\n");
	fprintf (init, "   vy = v_azimuth*(x-INCLINATION*z);\n");
	fprintf (init, "   vz = v_azimuth*INCLINATION*y;\n");
	fprintf (init, "   density = pow(radius,-SIGMASLOPE-1.)*interm1;\n");
	fprintf (init, "   energy = ASPECTRATIO*ASPECTRATIO/sqrt(radius*radius+1e-2);\n");
	fprintf (init, "   interm3 = interm1/(interm1+1e-3)*1.001;\n");
	fprintf (init, "   density = density*interm3 + (1.-interm3)*1e-2;\n");
	fprintf (init, "   energy  = energy *interm3 + (1.-interm3)*1.0;\n");
	break;
case cuve:
	fprintf (init, "   density = 1.0;\n");
	fprintf (init, "   energy = 1.0;\n");
	fprintf (init, "   x = radius*cos(azimuth);\n");
	fprintf (init, "   y = radius*sin(azimuth);\n");
	fprintf (init, "   interm1 = .5*(sqrt(3.)*(x-12.)+y);\n");
	fprintf (init, "   interm2 = .5*(-x+12. + sqrt(3.)*y);\n");
	fprintf (init, "   if ((fabs(interm1) < 3.) && (fabs(interm2) < 1.5)) {\n");
	fprintf (init, "   interm3 = .01*cos(interm1*2.*M_PI);\n");
	fprintf (init, "   density += interm3;\n");
	fprintf (init, "   vrad += -interm3*cos(azimuth-M_PI/3.);\n");
	fprintf (init, "   v_azimuth -= -interm3*sin(azimuth-M_PI/3.)/radius;\n");
	fprintf (init, "   }\n");
	fprintf (init, "   interm1 = .5*(sqrt(3.)*(x-12.)-y);\n");
	fprintf (init, "   interm2 = .5*(x-12. + sqrt(3.)*y);\n");
	fprintf (init, "   if ((fabs(interm1) < 3.) && (fabs(interm2) < 1.5)) {\n");
	fprintf (init, "   interm3 = .01*cos(interm1*2.*M_PI);\n");
	fprintf (init, "   density += interm3;\n");
	fprintf (init, "   vrad += -interm3*cos(azimuth+M_PI/3.);\n");
	fprintf (init, "   v_azimuth -= -interm3*sin(azimuth+M_PI/3.)/radius;\n");
	fprintf (init, "   }\n");
	break;
default:
	break;
}
}
fclose (init);
}
