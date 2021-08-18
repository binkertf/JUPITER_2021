/********************************/
/*                              */
/* This file is created         */
/* automatically during         */
/* compilation. Do not edit.    */
/* See perl script              */
/* "potparser.pl" for details   */
/*                              */
/********************************/

#include "pot.h"
#include "potcodes.h"
#include "jupiter.h"

void DumpPotCode (filename)
char *filename;
{
    FILE *pot;
    pot = prs_open (filename);
    if (!EXTERNALPOTENTIAL) {
	fprintf (pot, "No external potential is applied to the system.\n");
	fclose (pot);
	return;
    }
    switch (PotentialCode) {
case vkoffring2d:
	fprintf (pot, "y = radius * sin(colatitude); /* y is a work variable */\n");
	fprintf (pot, "z = radius*cos(colatitude);\n");
	fprintf (pot, "pot = -1.0/y+.5*1e-3/0.15*(z-1.0)*(z-1.0);\n");
	break;
case isoatm:
	fprintf (pot, "pot = z+1e5;\n");
	break;
case isoatm2d:
	fprintf (pot, "pot = 2.*x+1.*y+1e5;\n");
	break;
case planet3d:
	fprintf (pot, "pot = keplerian_pot(radius, azimuth, colatitude);\n");
	break;
case planet2d:
	fprintf (pot, "pot = keplerian_pot(radius, azimuth, colatitude);\n");
	break;
case rope:
	fprintf (pot, "pot = -1./radius;\n");
	fprintf (pot, "x = radius*cos(azimuth)*sin(colatitude);\n");
	fprintf (pot, "pot += -PLANETMASS \\n");
	fprintf (pot, "/ sqrt(1. + radius*radius -2.*radius*cos(azimuth) + SMOOTHING*SMOOTHING);\n");
	fprintf (pot, "pot += PLANETMASS*x;\n");
	break;
case kepler:
	fprintf (pot, "pot = -1.0/radius;\n");
	break;
case null:
	fprintf (pot, "pot = 0.0;\n");
	break;
case homogeneity1:
	fprintf (pot, "pot = -1.0/radius;\n");
	break;
case homogeneity2:
	fprintf (pot, "pot = -2.0e30*6.67e-11/radius;\n");
	break;
case GalacticCyl:
	fprintf (pot, "pot = 4.0e4*log(radius)+50./70./70.*z*z;\n");
	break;
case plan1d:
	fprintf (pot, "pot = SMOOTHING*log(z)+1e5;\n");
	break;
case star1d:
	fprintf (pot, "pot = SMOOTHING*.5*z*z+1e5;\n");
	break;
case planet2d_noind:
	fprintf (pot, "pot = -1./radius;\n");
	fprintf (pot, "pot -= PLANETMASS\\n");
	fprintf (pot, "/(sqrt(1.+radius*radius-2*radius*cos(azimuth)+SMOOTHING*SMOOTHING));\n");
	break;
case planetersatzt_norad:
	fprintf (pot, "pot = -1./radius;\n");
	fprintf (pot, "pot -= PLANETMASS\\n");
	fprintf (pot, "/(sqrt(azimuth*azimuth+SMOOTHING*SMOOTHING));\n");
	break;
case planet2d_trunc:
	fprintf (pot, "pot = -1./radius;\n");
	fprintf (pot, "x = sqrt(1.+radius*radius-2.*radius*cos(azimuth));\n");
	fprintf (pot, "if (x > 0.02) x = 0.02;\n");
	fprintf (pot, "pot -= PLANETMASS/sqrt(x*x+SMOOTHING*SMOOTHING);\n");
	break;
case planet_incli:
	fprintf (pot, "work1 = cos(t);\n");
	fprintf (pot, "yp = sin(t);\n");
	fprintf (pot, "xp = radius*cos(INCLINATION);\n");
	fprintf (pot, "zp = radius*sin(INCLINATION);\n");
	fprintf (pot, "x = radius*cos(azimuth)*sin(colatitude);\n");
	fprintf (pot, "y = radius*sin(azimuth)*sin(colatitude);\n");
	fprintf (pot, "z = radius*cos(colatitude);\n");
	fprintf (pot, "pot = -1./radius;\n");
	fprintf (pot, "pot -= PLANETMASS/sqrt((x-xp)*(x-xp)+(y-yp)*(y-yp)+(z-zp)*(z-zp)+SMOOTHING*SMOOTHING);\n");
	break;
case tide:
	fprintf (pot, "work1 = (0.01-OMEGAFRAME)*t;\n");
	fprintf (pot, "work2 = cos(azimuth-work1)*radius*sin(colatitude);\n");
	fprintf (pot, "pot = work2*work2*1e-5;\n");
	break;
case throop:
	fprintf (pot, "pot = -1.0/sqrt(x*x+y*y+z*z+1e-2);\n");
	break;
default:
	break;
}
fclose (pot);
}
