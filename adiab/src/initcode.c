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
addinitcode ("standingshock", standingshock);
addinitcode ("rotball", rotball);
addinitcode ("sod1diso", sod1diso);
addinitcode ("sod1d", sod1d);
addinitcode ("sod1dL", sod1dL);
addinitcode ("trivial", trivial);
addinitcode ("sod1da", sod1da);
addinitcode ("sod1dadiso", sod1dadiso);
addinitcode ("aw1d", aw1d);
addinitcode ("aw1da", aw1da);
addinitcode ("sphaw1d", sphaw1d);
addinitcode ("daw1d", daw1d);
addinitcode ("awsph1dlat", awsph1dlat);
addinitcode ("vkring1d", vkring1d);
addinitcode ("soundwave1i", soundwave1i);
addinitcode ("soundwave1a", soundwave1a);
addinitcode ("soundwave1adiso", soundwave1adiso);
addinitcode ("soundwave2", soundwave2);
addinitcode ("sodpm06", sodpm06);
addinitcode ("dustpm06", dustpm06);
addinitcode ("ring3d", ring3d);
addinitcode ("tiltawsph2d", tiltawsph2d);
addinitcode ("tiltawsph2da", tiltawsph2da);
addinitcode ("vkoffring2d", vkoffring2d);
addinitcode ("sod2d", sod2d);
addinitcode ("sod2dadiso", sod2dadiso);
addinitcode ("sod2dacyl", sod2dacyl);
addinitcode ("sod2da", sod2da);
addinitcode ("sod2db", sod2db);
addinitcode ("sod2daiso", sod2daiso);
addinitcode ("sod2dbiso", sod2dbiso);
addinitcode ("soundcyl2di", soundcyl2di);
addinitcode ("slushsph2di", slushsph2di);
addinitcode ("soundsph2di", soundsph2di);
addinitcode ("soundsph2da", soundsph2da);
addinitcode ("sodsph3di", sodsph3di);
addinitcode ("sodsph3da", sodsph3da);
addinitcode ("kepler3d", kepler3d);
addinitcode ("keplerdust", keplerdust);
addinitcode ("keplerdust_diff", keplerdust_diff);
addinitcode ("keplergas_diff", keplergas_diff);
addinitcode ("essai", essai);
addinitcode ("tides", tides);
addinitcode ("homogeneity1", homogeneity1);
addinitcode ("homogeneity2", homogeneity2);
addinitcode ("flatatm", flatatm);
addinitcode ("flatatma", flatatma);
addinitcode ("flatatmpl", flatatmpl);
addinitcode ("isoatm", isoatm);
addinitcode ("isoatma", isoatma);
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
addinitcode ("throop", throop);
addinitcode ("W19_1_g", W19_1_g);
addinitcode ("W19_1_d", W19_1_d);
addinitcode ("W19_1_d_05", W19_1_d_05);
addinitcode ("W19_1_d_2", W19_1_d_2);
addinitcode ("W19_1_d_3", W19_1_d_3);
addinitcode ("W19_1_d_4", W19_1_d_4);
}
