#include "jupiter.h"
#include "calc_gamma.h"


/** calculates Gamma for cell at [0,0,0] for all times during simulation*/
void TestGammaSingleCell() {
    FluidWork *fw; /* the current fluid patch (aka fluid work) we consider */
    real *density_ptr, *temperature_ptr; /*we want to read out the density and temperature of the fluid work */
    fw = CurrentFluidPatch; /* access current fluid patch (global variable) */
    density_ptr = fw->Density; /*read in density (1D array) */
    temperature_ptr = fw->Temperature; /*read in temperature (1D array) */
    real test_gamma;

    /* look at gamma value just for first cell
    use to calculate gamma, then print to log file*/
    test_gamma = Gamma1(*(temperature_ptr+0), *(density_ptr+0)); /* automatically points to first index*/
    pInfo ("### Value of gamma %g, at position (0,0,0) at date %.15g\n", test_gamma, GlobalDate);
    /* printf ("Value of gamma %.15g, at position (0,0,0) at date %.15g\n", test_gamma, GlobalDate); */

}


/* modeled pInfo and looping through dimensions on ExPot.c (lines 156 - 165) and predictor_adiab.c (lies 6 - 23) */
/** calculates Gamma for each cell in each FluidPatch for all times during simulation*/
void TestGammaFluidPatch(){
    long i, j, k, m; /* dimensions through which we loop (i,j,k) and counter m */
    FluidWork *fw; /* the current fluid patch (aka fluid work) we consider */
    real *density_ptr, *temperature_ptr; /*we want to read out the density and temperature of the fluid work */
    long stride[3], gncell[3]; /* variables needed to loop through all cells */
    real test_gamma = 0.;
    real x = 0.,y= 0.,z= 0.;
    real *center[3];

    fw = CurrentFluidPatch; /* access current fluid patch (global variable) */
    density_ptr = fw->Density; /*read in density (1D array) */
    temperature_ptr = fw->Temperature; /*read in temperature (1D array) */
    getgridsize (fw->desc, gncell, stride);

    for (i = 0; i < 3; i++)	/* 3, not NDIM */
        center[i] = fw->desc->Center[i];


         /* loop through all values (so all cells) in current fluid patch */
         for (k = 0; k < gncell[2]; k++) { /* ghosts are included */
             for (j = 0; j < gncell[1]; j++) {
                 for (i = 0; i < gncell[0]; i++) {
                     m = i * stride[0] + j * stride[1] + k * stride[2];

                     /*	want to be able to read off positions of cells*/
                     if (fw->desc->CommSource[m]) {
                         x = center[_RAD_][m];
                         y = center[_AZIM_][m];
                         z = center[_COLAT_][m];

                     /* for each cell, read out temp and density, then use to calculate gamma for this cell, then print to log file*/
                     test_gamma = Gamma1(*(temperature_ptr + m), *(density_ptr + m));
                     pInfo("### Value of gamma %g, at position x %g, y %g, z %g, at date %.15g\n", test_gamma, x, y, z,
                           GlobalDate);
                     /* printf("Value of gamma %.15g, at position x %g, y %g, z %g, at date %.15g\n", test_gamma, x, y, z,
                           GlobalDate); */
                 }
             }
         }
     }
}




/* modeled pInfo and looping through dimensions on ExPot.c (lines 156 - 165) and predictor_adiab.c (lies 6 - 23) */
/** prints all cell centers to logfile*/
void FindAllCellCenters(){
    long i, j, k, m; /* dimensions through which we loop (i,j,k) and counter m */
    FluidWork *fw; /* the current fluid patch (aka fluid work) we consider */
    long stride[3], gncell[3]; /* variables needed to loop through all cells */
    real x = 0.,y= 0.,z= 0.;
    real *center[3];

    fw = CurrentFluidPatch; /* access current fluid patch (global variable) */
    getgridsize (fw->desc, gncell, stride);

    for (i = 0; i < 3; i++)	/* 3, not NDIM */
        center[i] = fw->desc->Center[i];


    /* loop through all values (so all cells) in current fluid patch */
    for (k = 0; k < gncell[2]; k++) { /* ghosts are included */
        for (j = 0; j < gncell[1]; j++) {
            for (i = 0; i < gncell[0]; i++) {
                m = i * stride[0] + j * stride[1] + k * stride[2];

                /*	want to be able to read off positions of cells*/
                if (fw->desc->CommSource[m]) {
                    x = center[_RAD_][m];
                    y = center[_AZIM_][m];
                    z = center[_COLAT_][m];

                    /* pInfo("(x,y,z), %g, %g, %g\n", x, y, z); */
                    printf("(x,y,z), %g, %g, %g\n", x, y, z);
                }
            }
        }
    }
}



