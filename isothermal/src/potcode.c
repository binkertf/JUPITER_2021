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
static boolean PotCodesSet = NO;
void setpotcodes () {
if (PotCodesSet) return;
PotCodesSet = YES;
addpotcode ("vkoffring2d", vkoffring2d);
addpotcode ("isoatm", isoatm);
addpotcode ("isoatm2d", isoatm2d);
addpotcode ("planet3d", planet3d);
addpotcode ("planet2d", planet2d);
addpotcode ("rope", rope);
addpotcode ("kepler", kepler);
addpotcode ("null", null);
addpotcode ("homogeneity1", homogeneity1);
addpotcode ("homogeneity2", homogeneity2);
addpotcode ("GalacticCyl", GalacticCyl);
addpotcode ("plan1d", plan1d);
addpotcode ("star1d", star1d);
addpotcode ("planet2d_noind", planet2d_noind);
addpotcode ("planetersatzt_norad", planetersatzt_norad);
addpotcode ("planet2d_trunc", planet2d_trunc);
addpotcode ("planet_incli", planet_incli);
addpotcode ("tide", tide);
addpotcode ("throop", throop);
}
