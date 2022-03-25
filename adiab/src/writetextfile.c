#include "jupiter.h"
#include "calc_gamma.h"




/* modeled pInfo and looping through dimensions on ExPot.c (lines 156 - 165) and predictor_adiab.c (lies 6 - 23) */
/** calculates Gamma for each cell in each FluidPatch for all times during simulation*/
void CalcGammaFluidPatch(FILE *filehandle){
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

    fprintf(filehandle, "#\n");
    fprintf(filehandle, "D %g\n", GlobalDate); /* format: date */

    for (i = 0; i < 3; i++)	/* 3, not NDIM */
        center[i] = fw->desc->Center[i];


    /* loop through all values (so all cells) in current fluid patch, ghost cells NOT included */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
        for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
            for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
                m = i*stride[0]+j*stride[1]+k*stride[2];

                /*	want to be able to read off positions of cells*/
                if (fw->desc->CommSource[m]) {
                    x = center[_RAD_][m];
                    y = center[_AZIM_][m];
                    z = center[_COLAT_][m];

                    /* for each cell, read out temp and density, then use to calculate gamma for this cell, then print to log file*/
                    test_gamma = Gamma1(*(temperature_ptr + m), *(density_ptr + m));
                    fprintf(filehandle, "%g\n", test_gamma); /* format: GAMMA */
                }
            }
        }
    }
}




/* modeled pInfo and looping through dimensions on ExPot.c (lines 156 - 165) and predictor_adiab.c (lies 6 - 23) */
/** prints all cell centers to logfile*/

void FindAllCellCenters(FILE *filehandle){
    long i, j, k, m; /* dimensions through which we loop (i,j,k) and counter m */
    FluidWork *fw; /* the current fluid patch (aka fluid work) we consider */
    long stride[3], gncell[3]; /* variables needed to loop through all cells */
    real x = 0.,y= 0.,z= 0.;
    real *center[3];

    fw = CurrentFluidPatch; /* access current fluid patch (global variable) */
    getgridsize (fw->desc, gncell, stride);

    for (i = 0; i < 3; i++)	/* 3, not NDIM */
        center[i] = fw->desc->Center[i];


    /* loop through all values (so all cells) in current fluid patch, ghost cells NOT included */
    for (k = Nghost[2]; k < gncell[2]-Nghost[2]; k++) {
        for (j = Nghost[1]; j < gncell[1]-Nghost[1]; j++) {
            for (i = Nghost[0]; i < gncell[0]-Nghost[0]; i++) {
                m = i*stride[0]+j*stride[1]+k*stride[2];

                /*	want to be able to read off positions of cells*/
                if (fw->desc->CommSource[m]) {
                    x = center[_RAD_][m];
                    y = center[_AZIM_][m];
                    z = center[_COLAT_][m];

                    fprintf(filehandle,"%g, %g, %g\n", x,y, z);


                    /* pInfo("(x,y,z), %g, %g, %g\n", x, y, z); */
                    /* printf("(x,y,z), %g, %g, %g\n", x, y, z);*/

                }
            }
        }
    }
}


void WriteToFile(){
    if(GlobalDate==0) {
        fprintf(fpt, "%ld\n",SIZE1*SIZE2*SIZE3);
        fprintf(fpt, "CELL CENTERS\n");
        FindAllCellCenters(fpt);
    }
    CalcGammaFluidPatch(fpt);
}



