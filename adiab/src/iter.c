#include "jupiter.h"

real LevelDate[REFINEUPPERLIMIT];
extern long TimeStepRatio[REFINEUPPERLIMIT];
extern long BaseStepRatio[REFINEUPPERLIMIT];
extern boolean LevelHasChangedSinceCFL[REFINEUPPERLIMIT];

void ItereLevel (dt, level)
     real dt;
     long level;
{
    long printingCounter=0;
  tGrid_CPU *item;
  FluidPatch *Fluid;
  item = Grid_CPU_list;
  /* a TrueBC (somelev) must ALWAYS follow an ExecCommSame (somelev) */

  ExecCommSame (level-1); 	/* This has been tested with various */
  TrueBC (level-1);		/* grid configurations and CPU numbers */

  ExecCommUp (level-1);

  /* a TrueBC (somelev) must ALWAYS follow an ExecCommSame (somelev) */
  ExecCommSame (level);
  TrueBC (level);

  if (Stellar) {
    GetOpticalDepthFromLevel (level-1);
    GetEnergyRadFromLevel (level-1);
  }


  while (item != NULL) {
    if ((level == item->level) && (item->cpu == CPU_Rank)) {
      Fluid = item->Fluid;
      SetFluidProperties (Fluid);
      while (Fluid != NULL) {
          SendToCurrent(Fluid);

#if (TESTING_GAMMA==YES)
          if (GlobalDate == 0) {
              /* TestGammaValue(); */
          }

          /* !Isothermal*/
          if (!Isothermal) {
              double tchunk = DT * NTOT / 5.0;
              double cond1 = fabs(GlobalDate);
              double cond2 = fabs(GlobalDate - tchunk);
              double cond3 = fabs(GlobalDate - tchunk * 2);
              double cond4 = fabs(GlobalDate - tchunk * 3);
              double cond5 = fabs(GlobalDate - tchunk * 4);
              double epsilon = 1e-3;
              if (cond1 < epsilon || cond2 < epsilon || cond3 < epsilon || cond4 < epsilon || cond5 < epsilon) {
                  if (printingCounter == 0) {
                      /* prints only NTOT times, so once for every course timestep DT */
                      /*
                      TestGammaSingleCell();
                      TestGammaFluidPatch();
                      */
                      WriteToFile();
                      printingCounter++;
                  }
              }
          }

#endif



        if (Fluid->next==NULL){ //gas
          HydroKernel (dt);
        }else{ //dust
          if (diffmode != 2){
            DustKernel (dt);
          }else{
            SendToSecondary(Fluid->next);//create a second fluid patch for the gas to access in the diffusion calculation
            DustDiffKernel (dt);
          }
        }
        CurrentToPatch (Fluid);
	      Fluid = Fluid->next;
      }
    }
    FluidCoupling (item, dt);	/* Couple all fluids accessible to this item */
    MultifluidDiffusionPressure (item, dt); // Communicate turbulent diffusion pressure to dust fluids 
    /* This must be done after the individual steps on each fluid */
    item = item->next;
  }
    LevelHasChangedSinceCFL[level] = YES;
  ExecCommDownMean (level+1);
}

real RecursiveIteration (dt, level)
     real dt;
     long level;
{
  real time_var, dt_cfl_loc, dt_cfl;
  long i;
    if (level == 0) {		/* We begin a global timestep */
    for (i = 0; i <= LevMax; i++) {
      LevelHasChangedSinceCFL[i] = YES;
    }
    FindBestSubcycling ();
  }
  if (level < LevMax) { /* Not finest level */
    ResetFaceFluxesLevel (level+1);
    ExecCommDownMean (level+1);
    RecursiveIteration (dt/(real)TimeStepRatio[level], level+1);
      if (TimeStepRatio[level] > 1)
      RecursiveIteration (dt/(real)TimeStepRatio[level], level+1);
    ExecCommDownFlux (level+1);
    time_var = LevelDate[level+1]-LevelDate[level];
    ItereLevel (time_var, level);
      LevelDate[level] = LevelDate[level+1];
    return time_var;
  } else {			/* Finest level */
    dt_cfl_loc = CourantLimitGlobal () / (real)BaseStepRatio[level];
    MPI_Allreduce (&dt_cfl_loc, &dt_cfl, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    dt = (dt_cfl < dt ? dt_cfl : dt);
    ItereLevel (dt, level);
    LevelDate[level] += dt;
    return dt;
  }
}
