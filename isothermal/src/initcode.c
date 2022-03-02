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
static boolean InitCodesSet = NO;
void setinitcodes () {
if (InitCodesSet) return;
InitCodesSet = YES;
addinitcode ("multi1d", multi1d);
addinitcode ("multi1db", multi1db);
addinitcode ("multi1db_dust", multi1db_dust);
addinitcode ("W19_1_g", W19_1_g);
addinitcode ("W19_1_d", W19_1_d);
addinitcode ("W19_2_g", W19_2_g);
addinitcode ("W19_2_d", W19_2_d);
addinitcode ("sod1d", sod1d);
addinitcode ("sod1db", sod1db);
addinitcode ("aw1d", aw1d);
addinitcode ("sphaw1d", sphaw1d);
addinitcode ("daw1d", daw1d);
addinitcode ("awsph1dlat", awsph1dlat);
addinitcode ("vkring1d", vkring1d);
addinitcode ("soundwave1", soundwave1);
addinitcode ("soundwave2", soundwave2);
addinitcode ("sodpm06", sodpm06);
addinitcode ("dustpm06", dustpm06);
addinitcode ("tiltawsph2d", tiltawsph2d);
addinitcode ("vkoffring2d", vkoffring2d);
addinitcode ("sod2d", sod2d);
addinitcode ("kepler3d", kepler3d);
addinitcode ("keplerdust", keplerdust);
addinitcode ("keplerdust_2", keplerdust_2);
addinitcode ("essai", essai);
addinitcode ("tides", tides);
addinitcode ("homogeneity1", homogeneity1);
addinitcode ("homogeneity2", homogeneity2);
addinitcode ("flatatm", flatatm);
addinitcode ("flatatmpl", flatatmpl);
addinitcode ("isoatm", isoatm);
addinitcode ("isoatmpert", isoatmpert);
addinitcode ("isoatmpl", isoatmpl);
addinitcode ("isoatmpll", isoatmpll);
addinitcode ("isoatmps", isoatmps);
addinitcode ("isoatm2d", isoatm2d);
addinitcode ("isoatm2dp", isoatm2dp);
addinitcode ("trivialsphere", trivialsphere);
addinitcode ("HVCCyl", HVCCyl);
addinitcode ("testrot", testrot);
addinitcode ("testrefleq", testrefleq);
addinitcode ("plan1d", plan1d);
addinitcode ("star1d", star1d);
addinitcode ("local", local);
addinitcode ("localgrad", localgrad);
addinitcode ("novortgrad", novortgrad);
addinitcode ("sod2d_vert", sod2d_vert);
addinitcode ("throop", throop);
addinitcode ("cuve", cuve);
}
